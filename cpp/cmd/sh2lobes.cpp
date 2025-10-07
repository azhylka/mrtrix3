#include "command.h"
#include "dwi/fmls.h"
#include "command.h"
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

void usage ()
{
  AUTHOR = "Andrey Zhylka (zhylka.ai@gmail.com)";

  SYNOPSIS = "Extract the peaks' lobes of a spherical harmonic function in each voxel";

  DESCRIPTION
  + "Peaks of the spherical harmonic function in each voxel are located by "
    "commencing a Newton search along each of a set of pre-specified directions";

  DESCRIPTION
  + Math::SH::encoding_description;

  ARGUMENTS
  + Argument ("SH", "the input image of SH coefficients.")
    .type_image_in()
  + Argument ("fixel_directory", "the output fixel directory")
    .type_directory_out();


  OPTIONS
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
        H (header), fixel_count (0), max_per_voxel (maxnum), dir_from_peak (dir_from_peak) { }

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
            
            double lobe_width = 0;
            const Eigen::Vector3f& reference_peak(dir_from_peak ? lobe.get_peak_dir(0).cast<float>() : lobe.get_mean_dir().cast<float>());

            for (index_type peak_idx = 0; peak_idx != lobe.num_peaks(); peak_idx++) {
              auto lobe_dir = lobe.get_peak_dir(peak_idx).cast<float>();
              double deviation = atan2(reference_peak.cross(lobe_dir).norm(), reference_peak.dot(lobe_dir));
              lobe_width = max(lobe_width, 2*deviation);
            }
            
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
    index_type max_per_voxel;
    bool dir_from_peak;
};




bool Segmented_FOD_receiver::operator() (const DWI::FMLS::FOD_lobes& in)
{
  if (in.size()) {
    lobes.emplace_back (in, max_per_voxel, dir_from_peak);
    fixel_count += lobes.back().size();
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

  auto fixel_data_header (H);
  fixel_data_header.ndim() = 3;
  fixel_data_header.size(0) = fixel_count;
  fixel_data_header.size(2) = 1;
  fixel_data_header.datatype() = DataType::Float32;
  fixel_data_header.datatype().set_byte_order_native();

  if (dir_path.size()) {
    dir_image = make_unique<DataImage> (DataImage::create (Path::join(fixel_directory_path, dir_path), fixel_data_header));
  }
  if (peak_amp_path.size()) {
    peak_amp_image = make_unique<DataImage> (DataImage::create (Path::join(fixel_directory_path, peak_amp_path), fixel_data_header));
  }
  if (lobe_width_path.size()) {
    lobe_width_image = make_unique<DataImage> (DataImage::create (Path::join(fixel_directory_path, lobe_width_path), fixel_data_header));
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
      for (size_t i = 0; i < n_vox_fixels; ++i) {
        peak_amp_image->index(0) = offset + i;
        peak_amp_image->value() = vox_fixels[i].max_peak_amp;
      }
    }

    if (lobe_width_image) {
      for (size_t i = 0; i < n_vox_fixels; ++i) {
        lobe_width_image->index(0) = offset + i;
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
  if (!opt.empty())
    mask_data = Image<float>::open (std::string(opt[0][0]));
  
  DWI::Directions::FastLookupSet dirs (1281);
  
  DWI::FMLS::FODQueueWriter writer (fod, mask_data);
  DWI::FMLS::Segmenter fmls (dirs, Math::SH::LforN (fod.size(3)));
  fmls.set_integral_threshold(0.0);
  fmls.set_peak_value_threshold(0.0);
  
  const bool dir_as_peak = !get_options("dirpeak").empty();
  const index_type maxnum = get_option_value("maxnum", 0);

  Segmented_FOD_receiver receiver (H, maxnum, dir_as_peak);
  auto& fixel_directory_path  = argument[1];

  receiver.set_fixel_directory_output (fixel_directory_path);

  std::string file_extension (".mif");
  if (get_options ("nii").size())
    file_extension = ".nii";

  static const std::string default_index_filename ("index" + file_extension);
  static const std::string default_directions_filename ("directions" + file_extension);
  static const std::string lobe_width_filename ("lobe_width" + file_extension);
  static const std::string peak_amp_filename ("peak_amp" + file_extension);
  receiver.set_index_output(default_index_filename);
  receiver.set_directions_output(default_directions_filename);
  receiver.set_lobe_width_output(lobe_width_filename);
  receiver.set_peak_amp_output(peak_amp_filename);
  
  Thread::run_queue (writer, Thread::batch (DWI::FMLS::SH_coefs()), Thread::multi (fmls), Thread::batch (DWI::FMLS::FOD_lobes()), receiver);
  receiver.commit();
  
}

