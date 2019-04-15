layout(location = 0) in vec3 vertex;
layout(location = 1) in vec3 position;
layout(location = 2) in float vert_scale;

layout(std140, row_major) uniform viewMat{
    mat4 vpMatrix;
    mat4 rMatrix;
};
uniform mat3 pos_scale;
uniform float scale_fac;
uniform vec3 offset;
uniform vec4 color;

out vec3 normals_cameraspace;
out vec3 EyeDirection_cameraspace;
out vec3 LightDirection_cameraspace;
out vec4 MaterialDiffuseColor;

void main(void)
{
    gl_Position = vpMatrix * vec4(vertex * vert_scale * scale_fac +
                                  position * pos_scale +
                                  offset, 1);
    vec3 vertex_cameraspace = (rMatrix * vec4(vertex, 1)).xyz;
    EyeDirection_cameraspace = vec3(0,0,25) - vertex_cameraspace;
    LightDirection_cameraspace = vec3(10,10,10) + EyeDirection_cameraspace;
    normals_cameraspace = (rMatrix * vec4(vertex, 0)).xyz;
    MaterialDiffuseColor = color;
}

