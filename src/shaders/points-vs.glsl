#version 330
layout(location = 0) in vec3 a_Position;

uniform mat4 u_ModelViewProjectionMatrix;

void main() {
  gl_Position = u_ModelViewProjectionMatrix * vec4(a_Position, 1.0);
  gl_PointSize = 25.0 / gl_Position.w;
}
