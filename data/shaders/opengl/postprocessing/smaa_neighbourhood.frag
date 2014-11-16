
#ifndef SMAA_PIXEL_SIZE
#define SMAA_PIXEL_SIZE vec2(1.0 / 800.0, 1.0 / 600.0)
#endif
#define SMAA_PRESET_ULTRA 1
#define SMAA_GLSL_3 1
#define SMAA_ONLY_COMPILE_PS 1

#include "SMAA.glsl"

uniform sampler2D albedo_tex;
uniform sampler2D blend_tex;
in vec2 texcoord;
in vec4 offset[2];
in vec4 dummy2;

out vec4 frag_color;

void main()
{
  frag_color = SMAANeighborhoodBlendingPS(texcoord, offset, albedo_tex, blend_tex);
}