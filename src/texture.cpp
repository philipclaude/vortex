#include "texture.h"

#include <fmt/format.h>

#include <algorithm>

#include "log.h"

// TODO move these to header file
extern "C" {
unsigned char *stbi_load(char const *filename, int *x, int *y,
                         int *channels_in_file, int desired_channels);
void stbi_image_free(void *retval_from_stbi_load);
}

namespace vortex {

Texture::Texture(const std::string &filename, TextureOptions options)
    : channels_(kFormat2Channels.at(options.format)), flipy_(options.flipy) {
  read(filename);
}

void Texture::read(const std::string &filename) {
  // From the stb documentation (N is the number of desired channels):
  // An output image with N components has the following components interleaved
  // in this order in each pixel:
  //
  //     N=#comp     components
  //       1           grey
  //       2           grey, alpha
  //       3           red, green, blue
  //       4           red, green, blue, alpha
  int n;
  unsigned char *pixels;
  pixels = stbi_load(filename.c_str(), &width_, &height_, &n, channels_);
  ASSERT(pixels);

  // save the pixel data
  data_.resize(width_ * height_ * channels_);
  for (int i = 0; i < width_ * height_ * channels_; i++) data_[i] = pixels[i];
  stbi_image_free(pixels);

  LOG << fmt::format("read {} x {} image with {} channels (saved {})", width_,
                     height_, n, channels_);
}

void Texture::sample(double s, double t, double *value) const {
  ASSERT(s >= 0 && s <= 1);
  ASSERT(t >= 0 && t <= 1);
  int i = std::clamp(int(s * width_), 0, width_ - 1);
  int j = std::clamp(int(t * height_), 0, height_ - 1);
  if (flipy_) j = height_ - j - 1;
  for (int d = 0; d < channels_; d++) {
    value[d] = data_[(j * width_ + i) * channels_ + d];
  }
}

}  // namespace vortex