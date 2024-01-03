#version 330
layout (location = 0) in vec3 a_Position;

uniform mat4 u_ModelViewProjectionMatrix;
uniform mat4 u_ModelViewMatrix;
uniform mat4 u_NormalMatrix;

out vec3 v_Position;
out vec3 v_Normal;
noperspective out vec3 v_Point;

void main() {
  gl_Position = u_ModelViewProjectionMatrix * vec4(a_Position, 1.0);
  v_Position = (u_ModelViewMatrix * vec4(a_Position, 1.0)).xyz;
  v_Point = a_Position.xyz;
  v_Normal = mat3(u_NormalMatrix) * normalize(a_Position);
}
