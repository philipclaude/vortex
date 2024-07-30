#version 330

layout(location = 0) out vec4 fragColor;

in float v_Density;
in float v_DotProduct;

uniform samplerBuffer colormap;

#define ncolor 256

uniform float u_umin = 0.0;
uniform float u_umax = 1.0;

void main() {
  // Cull fragments where the dot product is less than or equal to zero
  if (v_DotProduct <= 0.0) {
    discard;
  }
  
  vec3 color;

  if (v_Density >= 0.0) {
    // Normalize the density value to the range [0, 1]
    float normalized_density = (v_Density - u_umin) / (u_umax - u_umin);
    normalized_density = clamp(normalized_density, 0.0, 1.0);

    // Map the normalized density to a colormap index
    int color_index = int(normalized_density * (ncolor - 1));

    // Fetch the color from the colormap
    color = texelFetch(colormap, color_index).rgb;
  } else {
    // Set color to black if no density information
    color = vec3(0.0, 0.0, 0.0);
  }

  float alpha = 1 / (1 + pow(5 * length(gl_PointCoord - vec2(0.5, 0.5)), 4));
  if (alpha < 0.5) discard;
  fragColor = vec4(color, alpha);
}
