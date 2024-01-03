#version 330

layout(location = 0) out vec4 fragColor;

uniform int u_lighting;
uniform int u_image;

uniform vec3 u_eye;
uniform float u_fov;
uniform int u_width;
uniform int u_height;

uniform mat4 u_InverseModel;
uniform mat4 u_Camera;

uniform sampler2D image;
uniform sampler2D normalmap;

// sphere parameters
vec3 center = vec3(0, 0, 0);
float radius = 1.0;

#define M_PI 3.141592654

vec3 get_color(in vec3 r, in vec3 eye) {
  // (eye - center) in object space
  vec3 U = eye - center;
  float B = dot(r, U);
  float C = length(U) * length(U) - radius * radius;
  float disc = B * B - C;
  if (disc < 0.0) discard;
  float t = -B - sqrt(disc);
  if (t < 0) discard;

  // surface point (p) and normal (n) in object space
  vec3 p = eye + t * r;
  vec3 n = normalize(p - center);

  // spherical and texture coordinates
  float theta = atan(p.y, p.x) + M_PI; // in [0, 2 * \pi]
  float phi = M_PI - acos(p.z); // in [0, \pi]
  float u = 0.5 * theta / M_PI;
  float v = phi / M_PI;

  // differentiating p = (x, y, z) to calculate dp/du
  // x = \cos(\theta) * sin(\phi)
  // y = \sin(\theta) * sin(\phi)
  // z = \cos(\phi);
  vec3 dp_du = vec3(-p.y, p.x, 0);
  vec3 dp_dv = cross(n, dp_du);

  vec3 km = vec3(0.5) * (1 - u_image) + u_image * texture(image, vec2(u, v)).rgb;

  // normal in tangent space
  vec3 nt = texture(normalmap, vec2(u, v)).xyz * 2.0 - 1.0;

  // map normal to object space
  n = normalize(nt.x * dp_du + nt.y * dp_dv + nt.z * n);

  vec3 ca = vec3(0.4); // ambient light color
  vec3 l = normalize(eye - p); // direction to light (eye) in object space
  return km * (ca + u_lighting * abs(dot(l, n)) + (1 - u_lighting) * 0.8);
}

void main() {
  // pixel coordinates relative to camera
  float height = 2.0 * tan(u_fov / 2.0);
  float width = height * (float(u_width) / float(u_height));
  float x = -0.5 * width + width * gl_FragCoord.x / u_width;
  float y = -0.5 * height + height * gl_FragCoord.y / u_height;

  // ray direction in object space
  vec3 r = mat3(u_InverseModel) * mat3(u_Camera) * vec3(x, y, -1);

  // intersect with sphere and determine color
  vec3 color = get_color(normalize(r), u_eye);
  fragColor = vec4(color, 1);
}
