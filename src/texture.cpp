#include "texture.h"

#include <fmt/format.h>

#include <algorithm>

#include "log.h"

// TODO move these to header file
extern "C" {
unsigned char *stbi_load(char const *filename, int *x, int *y,
                         int *channels_in_file, int desired_channels);
void stbi_image_free(void *retval_from_stbi_load);
int stbi_write_jpg(char const *filename, int x, int y, int comp,
                   const void *data, int quality);
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

void Texture::limit(double min, double max, bool reverse) {
  ASSERT(channels_ == 1);
  for (int j = 0; j < height_; j++) {
    for (int i = 0; i < width_; i++) {
      int k = j * width_ + i;
      if (data_[k] < min) data_[k] = min;
      if (data_[k] > max) data_[k] = max;
      if (reverse) data_[k] = 255 - data_[k];
    }
  }
}

void Texture::make_binary(double threshold, double min, double max) {
  ASSERT(channels_ == 1);
  for (int j = 0; j < height_; j++) {
    for (int i = 0; i < width_; i++) {
      int k = j * width_ + i;
      if (data_[k] < threshold)
        data_[k] = min;
      else
        data_[k] = max;
    }
  }
}

void Texture::make_periodic() {
  for (int j = 0; j < height_; j++) {
    int k0 = j * width_;
    int k1 = j * width_ + width_ - 1;
    data_[k0] = data_[k1];
  }
}

void Texture::smooth(int n_iter) {
  ASSERT(channels_ == 1);
  for (int iter = 0; iter < n_iter; iter++) {
    for (int j = 1; j < height_ - 1; j++) {
      for (int i = 1; i < width_ - 1; i++) {
        int k = j * width_ + i;
        int r = k + 1;
        int l = k - 1;
        int t = k + width_;
        int b = k - width_;
        // data_[k] = 0.25 * double(data_[r] + data_[l] + data_[t] + data_[b]);
        data_[k] = 0.125 * double(data_[r] + data_[l] + data_[t] + data_[b] +
                                  data_[t + 1] + data_[t - 1] + data_[b + 1] +
                                  data_[b - 1]);
      }
    }
    for (int j = 0; j < height_; j++) {
      int k = j * width_;
      data_[k] = data_[k + 1];
      data_[k + width_ - 1] = data_[k + width_ - 2];
    }
  }
}

void Texture::write(const std::string &filename) const {
  stbi_write_jpg(filename.c_str(), width_, height_, 1, data_.data(), 100);
}

}  // namespace vortex