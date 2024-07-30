#version 330
layout(location = 0) in vec3 a_Position;
layout(location = 1) in float a_Density;

out float v_Density;
out float v_DotProduct;

uniform mat4 u_ModelViewProjectionMatrix;
uniform mat4 u_ModelViewMatrix;
uniform int u_hasDensity;

void main() {
  vec3 normal = normalize(a_Position);
  vec3 transformedNormal = normalize((u_ModelViewMatrix * vec4(normal, 0.0)).xyz);
  vec3 viewDirection = vec3(0.0, 0.0, 1.0);

  float dotProduct = dot(transformedNormal, viewDirection);
  v_DotProduct = dotProduct;
   
  gl_Position = u_ModelViewProjectionMatrix * vec4(normalize(a_Position), 1.0);
  gl_PointSize = 8.0 / gl_Position.w;
  if (u_hasDensity == 1) {
    v_Density = a_Density;
  } else {
    v_Density = -1.0;
  }
}