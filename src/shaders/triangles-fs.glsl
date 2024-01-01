#version 330

layout(location = 0) out vec4 fragColor;

in vec3 v_Position;
in vec3 v_Normal;
noperspective in vec3 v_Point;

uniform sampler2D image;

uniform samplerBuffer colormap;
uniform samplerBuffer field;
uniform samplerBuffer texcoord;
uniform usamplerBuffer index;
uniform samplerBuffer points;
uniform usamplerBuffer visibility;
#ifdef SPLIT_PRIMITIVES
uniform usamplerBuffer primitive2cell;
#endif

uniform vec3 constant_color;
uniform int use_constant_color;
uniform int u_lighting;
uniform int u_clip;
uniform vec3 u_clip_center;
uniform vec3 u_clip_normal;
uniform int u_field_mode;
uniform int u_picking = -1;

uniform mat4 u_ModelViewProjectionMatrix;
uniform mat4 u_ModelViewMatrix;
uniform mat4 u_NormalMatrix;
uniform vec2 u_ViewportSize = vec2(800, 600);

#define LARGE_DISTANCE 1000000

// TODO: make these uniforms
const int ncolor = 256;
uniform float u_umin = 0;
uniform float u_umax = 1;

uniform int u_edges;
uniform float u_alpha;

void get_color(float u, inout vec3 color, in int alpha) {
  float umin = u_umin;
  float umax = u_umax;
  int indx = clamp(int(ncolor * (u - umin) / (umax - umin)), 0, 255);
  color = (1 - alpha) * color + alpha * texelFetch(colormap, indx).xyz;
}

void shading(in vec3 l, in vec3 n, in vec3 color, out vec3 color_out) {
  vec3 cd = color * max(0.0, dot(l, n));
  vec3 ca = vec3(0.2) * color;
  color_out = ca + cd;
}

vec3 barycentric(in vec3 a, in vec3 b, in vec3 c, in vec3 p) {
  float ia = 1.0 / length(cross(b - a, c - a));
  float a0 = length(cross(c - b, p - b));
  float a1 = length(cross(p - a, c - a));
  float a2 = length(cross(b - a, p - a));
  return vec3(a0 * ia, a1 * ia, a2 * ia);
}

void main() {


  int id = gl_PrimitiveID; // triangle
#ifdef SPLIT_PRIMITIVES
  uint cell = texelFetch(primitive2cell, id).r;
#else
  int cell = id;
#endif
  uint t0 = texelFetch(index, 3 * id).r;
  uint t1 = texelFetch(index, 3 * id + 1).r;
  uint t2 = texelFetch(index, 3 * id + 2).r;

  vec3 x0 = texelFetch(points, int(t0)).xyz;
  vec3 x1 = texelFetch(points, int(t1)).xyz;
  vec3 x2 = texelFetch(points, int(t2)).xyz;

  vec3 uvw = barycentric(x0, x1, x2, v_Point);
  vec4 p0 = u_ModelViewProjectionMatrix * vec4(x0, 1.0);
  vec4 p1 = u_ModelViewProjectionMatrix * vec4(x1, 1.0);
  vec4 p2 = u_ModelViewProjectionMatrix * vec4(x2, 1.0);

  vec2 q0 = p0.xy / p0.w;
  vec2 q1 = p1.xy / p1.w;
  vec2 q2 = p2.xy / p2.w;

  vec2 v1 = u_ViewportSize * (q1 - q0);
  vec2 v2 = u_ViewportSize * (q2 - q0);
  vec2 v3 = u_ViewportSize * (q2 - q1);

  int vis1 = int(texelFetch(visibility, 3 * id).r);
  int vis2 = int(texelFetch(visibility, 3 * id + 1).r);
  int vis3 = int(texelFetch(visibility, 3 * id + 2).r);

  float area = abs(v1.x * v2.y - v1.y * v2.x); // twice the area
  float h1 = (1 + LARGE_DISTANCE * vis1) * area / length(v3);
  float h2 = (1 + LARGE_DISTANCE * vis2) * area / length(v2);
  float h3 = (1 + LARGE_DISTANCE * vis3) * area / length(v1);

  float a1 = h1 * uvw.x;
  float a2 = h2 * uvw.y;
  float a3 = h3 * uvw.z;
  float d = min(min(a1, a2), a3);

  float alpha = 1.0;
  float intensity = u_edges * exp2(-0.25 * d * d);

  vec3 position = v_Position;
  vec3 normal = normalize(v_Normal);

  vec3 color = vec3(0.8, 0.8, 0.8);

#if ORDER == 0

  float u = texelFetch(field, int(cell)).r;
  get_color(u, color, u_field_mode);

#elif ORDER == 1

  float u0 = texelFetch(field, 3 * id).r;
  float u1 = texelFetch(field, 3 * id + 1).r;
  float u2 = texelFetch(field, 3 * id + 2).r;

  float s = uvw.x;
  float t = uvw.y;

  float u = (1.0 - s - t) * u0 + s * u1 + t * u2;
  get_color(u, color, u_field_mode);

#endif

  vec3 color_shaded;
  vec3 light_dir = -normalize(position); // light at eye
  shading(light_dir, normal, color, color_shaded);

  vec3 color_out = color * (1 - u_lighting) + u_lighting * color_shaded;

  fragColor = intensity * vec4(0, 0, 0, 1.0) +
              (1.0 - intensity) * vec4(color_out, alpha);

  // modify the color if this element was picked
  // u_picking = -1 if not picking, u_picking = element id if picking
  if (u_picking == int(cell))
    fragColor = intensity * vec4(0, 0, 0, 1.0) +
                (1.0 - intensity) * vec4(0.2, 0.5, 1.0, 1.0);
}
