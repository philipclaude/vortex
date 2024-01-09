#version 330 core

layout (points) in;
uniform mat4 u_ModelViewProjectionMatrix;
uniform mat4 u_NormalMatrix;
uniform mat4 u_ModelViewMatrix;

uniform int u_clip;
uniform vec3 u_clip_point;
uniform vec3 u_clip_normal;
uniform int u_earth;

uniform vec2 u_aa_radius;
uniform int u_width;
uniform int u_height;

flat out int id;

layout (triangle_strip, max_vertices = 6) out;

float t = 0.025;

uniform samplerBuffer points;
uniform usamplerBuffer index;

flat in int v_id[];

void main() {

  vec2 ViewportSize = vec2(u_width, u_height);

  uint i0 = texelFetch(index, 2 * v_id[0]    ).r;
  uint i1 = texelFetch(index, 2 * v_id[0] + 1).r;

  vec3 x0 = texelFetch(points, int(i0)).xyz;
  vec3 x1 = texelFetch(points, int(i1)).xyz;

  vec3 x0v = mat3(u_NormalMatrix) * x0;
  if (u_earth > 0 && x0v.z < 0) return;

  vec4 p0 = u_ModelViewProjectionMatrix * vec4(x0, 1.0);
  vec4 p1 = u_ModelViewProjectionMatrix * vec4(x1, 1.0);

  float u_width = ViewportSize[0];
  float u_height = ViewportSize[1];
  float u_aspect_ratio = u_height / u_width;

  vec2 v_line_width;
  v_line_width[0] = 5;
  v_line_width[1] = 5;

  vec2 ndc_a = p0.xy / p0.w;
  vec2 ndc_b = p1.xy / p1.w;

  vec2 line_vector = ndc_b - ndc_a;
  vec2 viewport_line_vector = line_vector * ViewportSize;
  vec2 dir = normalize(vec2(line_vector.x, line_vector.y * u_aspect_ratio));
  
  float line_width_a = max(1.0, v_line_width[0]) + u_aa_radius[0];
  float line_width_b = max(1.0, v_line_width[1]) + u_aa_radius[0];
  float extension_length = u_aa_radius[1];
  float line_length = length( viewport_line_vector ) + 2.0 * extension_length;
  
  vec2 normal    = vec2( -dir.y, dir.x );
  vec2 normal_a  = vec2( line_width_a/u_width, line_width_a/u_height ) * normal;
  vec2 normal_b  = vec2( line_width_b/u_width, line_width_b/u_height ) * normal;
  vec2 extension = vec2( extension_length / u_width, extension_length / u_height ) * dir;

  gl_Position = vec4( (ndc_a + normal_a - extension) * p0.w, p0.zw );
  id = v_id[0];
  EmitVertex();
  
  gl_Position = vec4( (ndc_a - normal_a - extension) * p0.w, p0.zw );
  id = v_id[0];
  EmitVertex();

  gl_Position = vec4( (ndc_b + normal_b + extension) * p1.w, p1.zw );
  id = v_id[0];
  EmitVertex();
  
  gl_Position = vec4( (ndc_b - normal_b + extension) * p1.w, p1.zw );
  id = v_id[0];
  EmitVertex();
  
  EndPrimitive();
}
