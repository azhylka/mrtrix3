#include "command.h"
#include "dwi/fmls.h"
#include "command.h"
#include "algo/loop.h"
#include "dwi/directions/set.h"
#include "file/path.h"
#include "fixel/fixel.h"
#include "fixel/helpers.h"
#include "image.h"
#include "math/SH.h"
#include "progressbar.h"
#include "thread_queue.h"
#include "memory.h"


#include <cmath>
#include <vector>
#include <iostream>

using namespace MR;
using namespace App;
using namespace std;

using index_type = unsigned int;

// Fraction of a lobe's maximum amplitude below which a direction takes no part in the width
//   calculation; see Primitive_FOD_lobes::Primitive_FOD_lobes() for what this measures.
constexpr double lobe_width_amplitude_fraction = 0.25;

void usage ()
{
  AUTHOR = "Andrey Zhylka (zhylka.ai@gmail.com)";

  SYNOPSIS = "Extract the peaks' lobes of a spherical harmonic function in each voxel";

  DESCRIPTION
  + "Peaks of the spherical harmonic function in each voxel are located by "
    "commencing a Newton search along each of a set of pre-specified directions"

  + "The fixel directions themselves are written to the nominated fixel directory. "
    "Each of the scalar metrics, however, is written as a 4D image with the same spatial "
    "dimensions and transform as the input, the fourth axis indexing the fixels within each "
    "voxel; the fixels of a voxel are ordered by decreasing lobe integral, so that volume 0 "
    "corresponds to the largest lobe. The length of that axis is the largest number of fixels "
    "found in any voxel, and those voxels that contain fewer fixels are padded with NaN.";

  DESCRIPTION
  + Math::SH::encoding_description;

  ARGUMENTS
  + Argument ("SH", "the input image of SH coefficients.")
    .type_image_in()
  + Argument ("fixel_directory", "the output fixel directory")
    .type_directory_out();


  OPTIONS
  + OptionGroup ("Metric values for voxel-wise output images")

  + Option ("lobe_width",
            "output the angular width of the FOD lobe per fixel, in degrees"
            " (twice the largest angle between the fixel direction and any direction of the lobe"
            " at which the FOD amplitude is at least one quarter of the lobe's maximum),"
            " as a 4D image whose fourth axis indexes the fixels within each voxel")
  + Argument ("image").type_image_out()

  + Option ("peak_amp",
            "output the amplitude of the FOD at the maximal peak per fixel,"
            " as a 4D image whose fourth axis indexes the fixels within each voxel")
  + Argument ("image").type_image_out()

  + DWI::FMLS::FMLSSegmentOption

  + OptionGroup ("Other options for sh2lobes")

   + Option ("mask",
            "only perform computation within the specified binary brain mask image.")
  + Argument ("image").type_image_in()
  
  + Option ("maxnum", "maximum number of fixels to output for any particular voxel (default: no limit)")
  + Argument ("number").type_integer(1)

  + Option ("nii", "output the directions and index file in nii format (instead of the default mif)")

  + Option ("dirpeak", "define the fixel direction as that of the lobe's maximal peak as opposed to its weighted mean direction (the default)");

}


class Segmented_FOD_receiver { 

  public:
    Segmented_FOD_receiver (const Header& header, const index_type maxnum = 0, bool dir_from_peak = false) :
        H (header), fixel_count (0), max_fixels_per_voxel (0), max_per_voxel (maxnum), dir_from_peak (dir_from_peak) { }

    void commit ();

    void set_fixel_directory_output (const std::string& path) { fixel_directory_path = path; }
    void set_index_output (const std::string& path) { index_path = path; }
    void set_directions_output (const std::string& path) { dir_path = path; }
    void set_lobe_width_output (const std::string& path) { lobe_width_path = path; }
    void set_peak_amp_output (const std::string& path) { peak_amp_path = path; }

    bool operator() (const DWI::FMLS::FOD_lobes&);


  private:

    struct Primitive_FOD_lobe { 
      Eigen::Vector3f dir;
      float lobe_width;
      float max_peak_amp;
      Primitive_FOD_lobe (Eigen::Vector3f dir, float max_peak_amp, float lobe_width) :
          dir (dir), max_peak_amp (max_peak_amp), lobe_width (lobe_width) {}
    };


    class Primitive_FOD_lobes : public vector<Primitive_FOD_lobe> {
      public:
        Primitive_FOD_lobes (const DWI::FMLS::FOD_lobes& in, const index_type maxcount, bool dir_from_peak) :
            vox (in.vox)
        {
          const index_type N = maxcount ? std::min (index_type(in.size()), maxcount) : in.size();
          for (index_type i = 0; i != N; ++i) {
            const DWI::FMLS::FOD_lobe& lobe (in[i]);
            
            // The width of the lobe is the angular extent of its core about the fixel direction:
            //   twice the largest angle subtended between that direction and any direction of the
            //   lobe at which the FOD amplitude reaches a quarter of the lobe's maximum. Taking
            //   the lobe's full extent instead would measure the position of its zero crossing,
            //   which sits on the shallow tail of the profile and so responds more to the lobe's
            //   overall scale than to the shape of the peak.
            // The reference amplitude is the largest amplitude among the directions the segmenter
            //   assigned to the lobe, rather than the Newton-refined peak value reported by
            //   -peak_amp: both the threshold and the values tested against it are then drawn from
            //   the same sampling of the FOD, and the direction attaining the maximum is always
            //   retained, so a lobe always has a well-defined width.
            // Note that the lobe's directions are not confined to a single hemisphere, so it is
            //   the acute angle that must be measured in each case.
            const Eigen::Vector3d reference_dir (dir_from_peak ? lobe.get_peak_dir(0) : lobe.get_mean_dir());

            const DWI::Directions::Mask& lobe_mask (lobe.get_mask());
            const DWI::Directions::Set& dirs (lobe_mask.get_dirs());
            const Eigen::Array<default_type, Eigen::Dynamic, 1>& values (lobe.get_values());

            double max_amplitude = 0.0;
            for (size_t dir_idx = 0; dir_idx != lobe_mask.size(); ++dir_idx) {
              if (lobe_mask[dir_idx])
                max_amplitude = std::max (max_amplitude, double (values[dir_idx]));
            }
            const double amplitude_threshold = lobe_width_amplitude_fraction * max_amplitude;

            double lobe_width = 0;
            for (size_t dir_idx = 0; dir_idx != lobe_mask.size(); ++dir_idx) {
              if (!lobe_mask[dir_idx] || values[dir_idx] < amplitude_threshold)
                continue;
              const Eigen::Vector3d& lobe_dir (dirs[dir_idx]);
              const double deviation = std::atan2 (reference_dir.cross (lobe_dir).norm(), std::abs (reference_dir.dot (lobe_dir)));
              lobe_width = std::max (lobe_width, 2.0 * deviation);
            }
            lobe_width *= 180.0 / Math::pi;

            if (dir_from_peak)
              this->emplace_back (lobe.get_peak_dir(0).cast<float>(), lobe.get_max_peak_value(), lobe_width);
            else
              this->emplace_back (lobe.get_mean_dir().cast<float>(), lobe.get_max_peak_value(), lobe_width);
          }
        }
        Eigen::Array3i vox;
    };

    Header H;
    std::string fixel_directory_path, index_path, peak_amp_path, dir_path, lobe_width_path;
    vector<Primitive_FOD_lobes> lobes;
    index_type fixel_count;
    index_type max_fixels_per_voxel;
    index_type max_per_voxel;
    bool dir_from_peak;
};




bool Segmented_FOD_receiver::operator() (const DWI::FMLS::FOD_lobes& in)
{
  if (in.size()) {
    lobes.emplace_back (in, max_per_voxel, dir_from_peak);
    fixel_count += lobes.back().size();
    // Note that this counts the lobes retained after any -maxnum truncation, which is applied
    //   by the Primitive_FOD_lobes constructor; the option therefore bounds the length of the
    //   fourth axis of the metric images without any further handling.
    max_fixels_per_voxel = std::max (max_fixels_per_voxel, index_type (lobes.back().size()));
  }
  return true;
}



void Segmented_FOD_receiver::commit ()
{
  if (!lobes.size() || !fixel_count)
    return;

  using DataImage = Image<float>;
  using IndexImage = Image<index_type>;

  const auto index_filepath = Path::join (fixel_directory_path, index_path);

  std::unique_ptr<IndexImage> index_image;
  std::unique_ptr<DataImage> dir_image;
  std::unique_ptr<DataImage> peak_amp_image;
  std::unique_ptr<DataImage> lobe_width_image;

  auto index_header (H);
  index_header.keyval()[Fixel::n_fixels_key] = str(fixel_count);
  index_header.ndim() = 4;
  index_header.size(3) = 2;
  index_header.datatype() = DataType::from<index_type>();
  index_header.datatype().set_byte_order_native();
  index_image = make_unique<IndexImage> (IndexImage::create (index_filepath, index_header));

  if (dir_path.size()) {
    auto dir_header (H);
    dir_header.ndim() = 3;
    dir_header.size(0) = fixel_count;
    dir_header.size(1) = 3;
    dir_header.size(2) = 1;
    dir_header.transform().setIdentity();
    dir_header.spacing(0) = dir_header.spacing(1) = dir_header.spacing(2) = 1.0;
    dir_header.datatype() = DataType::Float32;
    dir_header.datatype().set_byte_order_native();
    dir_image = make_unique<DataImage> (DataImage::create (Path::join(fixel_directory_path, dir_path), dir_header));
  }

  // The scalar metrics are written as standalone 4D images that retain the spatial dimensions
  //   and transform of the input, the fourth axis indexing the fixels within each voxel. The
  //   length of that axis is the largest fixel count encountered; every voxel is initialised to
  //   NaN so that those holding fewer fixels, and those holding none at all, are padded rather
  //   than reporting a width or an amplitude of zero.
  auto metric_header (H);
  metric_header.ndim() = 4;
  metric_header.size(3) = max_fixels_per_voxel;
  metric_header.datatype() = DataType::Float32;
  metric_header.datatype().set_byte_order_native();

  if (peak_amp_path.size())
    peak_amp_image = make_unique<DataImage> (DataImage::create (peak_amp_path, metric_header));
  if (lobe_width_path.size())
    lobe_width_image = make_unique<DataImage> (DataImage::create (lobe_width_path, metric_header));

  for (DataImage* metric_image : {peak_amp_image.get(), lobe_width_image.get()}) {
    if (!metric_image)
      continue;
    for (auto l = Loop(*metric_image)(*metric_image); l; ++l)
      metric_image->value() = NaN;
  }

  size_t offset (0);
  for (const auto& vox_fixels : lobes) {
    size_t n_vox_fixels = vox_fixels.size();

    assign_pos_of (vox_fixels.vox).to (*index_image);

    index_image->index(3) = 0;
    index_image->value () = n_vox_fixels;

    index_image->index(3) = 1;
    index_image->value() = offset;

    if (dir_image) {
      for (size_t i = 0; i < n_vox_fixels; ++i) {
        dir_image->index(0) = offset + i;
        dir_image->row(1) = vox_fixels[i].dir;
      }
    }

    if (peak_amp_image) {
      assign_pos_of (vox_fixels.vox).to (*peak_amp_image);
      for (size_t i = 0; i < n_vox_fixels; ++i) {
        peak_amp_image->index(3) = i;
        peak_amp_image->value() = vox_fixels[i].max_peak_amp;
      }
    }

    if (lobe_width_image) {
      assign_pos_of (vox_fixels.vox).to (*lobe_width_image);
      for (size_t i = 0; i < n_vox_fixels; ++i) {
        lobe_width_image->index(3) = i;
        lobe_width_image->value() = vox_fixels[i].lobe_width;
      }
    }

    offset += n_vox_fixels;
  }

  assert (offset == fixel_count);
}



void run ()
{
  Header H = Header::open (argument[0]);
  Math::SH::check (H);
  auto fod = H.get_image<float>();

  check_3D_nonunity (fod);

  auto opt = get_options ("mask");
  Image<float> mask_data;
  if (!opt.empty()) {
    mask_data = Image<float>::open (std::string(opt[0][0]));
    if (!dimensions_match (fod, mask_data, 0, 3))
      throw Exception ("Cannot use image \"" + str(opt[0][0]) + "\" as mask image;"
                       " dimensions do not match FOD image");
  }

  DWI::Directions::FastLookupSet dirs (1281);

  DWI::FMLS::FODQueueWriter writer (fod, mask_data);
  DWI::FMLS::Segmenter fmls (dirs, Math::SH::LforN (fod.size(3)));
  DWI::FMLS::load_fmls_thresholds (fmls);

  const bool dir_as_peak = !get_options("dirpeak").empty();
  const index_type maxnum = get_option_value("maxnum", 0);

  Segmented_FOD_receiver receiver (H, maxnum, dir_as_peak);
  auto& fixel_directory_path  = argument[1];

  Fixel::check_fixel_directory (fixel_directory_path, true, true);

  receiver.set_fixel_directory_output (fixel_directory_path);

  std::string file_extension (".mif");
  if (get_options ("nii").size())
    file_extension = ".nii";

  static const std::string default_index_filename ("index" + file_extension);
  static const std::string default_directions_filename ("directions" + file_extension);
  receiver.set_index_output(default_index_filename);
  receiver.set_directions_output(default_directions_filename);

  opt = get_options ("lobe_width");
  if (!opt.empty())
    receiver.set_lobe_width_output (opt[0][0]);
  opt = get_options ("peak_amp");
  if (!opt.empty())
    receiver.set_peak_amp_output (opt[0][0]);

  Thread::run_queue (writer, Thread::batch (DWI::FMLS::SH_coefs()), Thread::multi (fmls), Thread::batch (DWI::FMLS::FOD_lobes()), receiver);
  receiver.commit();
  
}

