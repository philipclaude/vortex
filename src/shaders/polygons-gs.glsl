#version 330 core

layout (points) in;

uniform mat4 u_ModelViewProjectionMatrix;
uniform mat4 u_NormalMatrix;
uniform mat4 u_ModelViewMatrix;
uniform mat4 u_ViewportMatrix;

uniform float u_alpha;

uniform vec2 u_ViewportSize = vec2(800, 600);

noperspective out vec3 v_Position;
noperspective out vec3 altitude;

uniform samplerBuffer coordinates;
uniform usamplerBuffer index;
uniform usamplerBuffer first;
uniform usamplerBuffer count;

layout (triangle_strip , max_vertices = 50) out;

flat in int[] v_id;
flat out int id;

#define LARGE_DISTANCE 1000000 // anything larger than max(width, height) would be fine

void main() {

  int n = int(texelFetch(count, v_id[0]).r);
  int f = int(texelFetch(first, v_id[0]).r);

  // the following can be used for convex polygons, but the second version is better for general polygons
  uint i0 = texelFetch(index, f).r;
  uint i1 = texelFetch(index, f + 1).r;

  vec3 X0 = texelFetch(coordinates, int(i0) ).xyz;
  vec3 X1 = texelFetch(coordinates, int(i1) ).xyz;

  vec4 x0 = vec4(X0, 1.0);
  vec4 x1 = vec4(X1, 1.0);

  vec4 p0 = u_ModelViewProjectionMatrix * x0;
  vec4 p1 = u_ModelViewProjectionMatrix * x1;

  for (int i = 2; i < n; i++) {

    uint i2 = texelFetch(index, int(f + i)).r;
    vec4 x2 = vec4(texelFetch(coordinates, int(i2)).xyz, 1.0);
    vec4 p2 = u_ModelViewProjectionMatrix * x2;

    //vec4 u = x1 - x0;
    //vec4 v = x2 - x0;
    //vec3 normal = normalize(mat3(u_NormalMatrix) * cross(u.xyz,v.xyz));
    //if (floor(u_alpha) * normal[2] < 0) return; // normal is facing away, this will be false if there is any transparency active

    #if 1
    vec2 q0 = p0.xy / p0.w;
    vec2 q1 = p1.xy / p1.w;
    vec2 q2 = p2.xy / p2.w;

    vec2 v1 = u_ViewportSize * (q1 - q0);
    vec2 v2 = u_ViewportSize * (q2 - q0);
    vec2 v3 = u_ViewportSize * (q2 - q1);

    float area = abs(v1.x * v2.y - v1.y * v2.x); // twice the area
    #else
    vec3 q0 = (u_ModelViewMatrix * x0).xyz;
    vec3 q1 = (u_ModelViewMatrix * x1).xyz;
    vec3 q2 = (u_ModelViewMatrix * x2).xyz;

    vec3 v1 = q1 - q0;
    vec3 v2 = q2 - q0;
    vec3 v3 = q2 - q1;

    float area = cross(v1, v2).z;

    #endif

    float h1 = area / length(v3);
    float h2 = area / length(v2);
    float h3 = area / length(v1);

    // these are interior, so set the distance far away
    #if 0
    int beta = (i == 2 && n != 3) ? 1 : 0;
    h2 = beta * LARGE_DISTANCE + (1 - beta) * h2;
    beta = (i == n - 1 && n != 3) ? 1 : 0;
    h3 = beta * LARGE_DISTANCE + (1 - beta) * h3;
    //h2 = LARGE_DISTANCE;
    //h3 = LARGE_DISTANCE;
    #else
    if (i == 2 && n != 3) {
      h2 = LARGE_DISTANCE;
    }
    else if (i == n-1 && n != 3) {
      h3 = LARGE_DISTANCE;
    }
    else {
      h2 = LARGE_DISTANCE;
      h3 = LARGE_DISTANCE;
    }
    #endif

    gl_Position = p0;
    v_Position  = x0.xyz;
    altitude    = vec3(h1, 0, 0);
    id = v_id[0];
    EmitVertex();

    gl_Position = p1;
    v_Position  = x1.xyz;
    altitude    = vec3(0, h2, 0);
    id = v_id[0];
    EmitVertex();

    gl_Position = p2;
    v_Position  = x2.xyz;
    altitude    = vec3(0, 0, h3);
    id = v_id[0];

    gl_PrimitiveID = gl_PrimitiveIDIn;
    EmitVertex();

    EndPrimitive();

    p1 = p2;
    x1 = x2;
  }

}
