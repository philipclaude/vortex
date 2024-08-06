#version 330
layout(location = 0) in vec3 a_Position;
layout(location = 1) in float a_Density;

out float v_Density;
out vec3 v_Normal;
out vec3 v_Position;

uniform mat4 u_ModelViewProjectionMatrix;
uniform mat4 u_ModelViewMatrix;
uniform int u_hasDensity;
uniform mat4 u_InverseModel;
uniform mat4 u_Camera;
uniform float u_fov;
uniform int u_width;
uniform int u_height;
uniform vec3 u_eye;

// sphere parameters
vec3 center = vec3(0, 0, 0);
float radius = 1.0;

void main() {
  vec4 position = vec4(a_Position, 1.0);

  // Calculate normal in view space
  vec3 normal = normalize((u_ModelViewMatrix * vec4(a_Position, 0.0)).xyz);
  v_Normal = normal;

  v_Position = a_Position;

  // Transform the position to clip space
  vec4 clip_position = u_ModelViewProjectionMatrix * position;
  vec4 view_position = u_ModelViewMatrix * position;
  
  // Convert from clip space to NDC
  vec3 ndc = clip_position.xyz / clip_position.w;

  // Calculate pixel coordinates in NDC space
  float height = 2.0 * tan(u_fov / 2.0);
  float width = height * (float(u_width) / float(u_height));
  float x = ndc.x * width / 2.0;
  float y = ndc.y * height / 2.0;

  // ray direction in object space
  vec3 r = normalize(mat3(u_InverseModel) * mat3(u_Camera) * vec3(x, y, -1));

  // (eye - center) in object space
  vec3 U = u_eye - center;
  float B = dot(r, U);
  float C = length(U) * length(U) - radius * radius;
  float disc = B * B - C;
  
  if (disc < 0.0) {
      position = vec4(2.0, 2.0, 2.0, 1.0); 
  } else {
      float t = -B - sqrt(disc);
      
      // If both intersections are behind the ray origin, move the point out of view
      if (t < 0.0) {
          position = vec4(2.0, 2.0, 2.0, 1.0); 
      } else {
          vec3 p = u_eye + t * r;
          if (length(p - a_Position) > 0.1) {
              position = vec4(2.0, 2.0, 2.0, 1.0);
          }
      }
  }
  gl_Position = u_ModelViewProjectionMatrix * position;
  gl_PointSize = 10.0 / gl_Position.w;

  if (u_hasDensity == 1) {
    v_Density = a_Density;
  } else {
    v_Density = -1.0;
  }
}