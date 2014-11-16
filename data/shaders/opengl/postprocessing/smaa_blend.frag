
#ifndef SMAA_PIXEL_SIZE
#define SMAA_PIXEL_SIZE vec2(1.0 / 800.0, 1.0 / 600.0)
#endif
#define SMAA_PRESET_ULTRA 1
#define SMAA_GLSL_3 1
#define SMAA_ONLY_COMPILE_PS 1

#include "SMAA.glsl"

uniform sampler2D edge_tex;
uniform sampler2D area_tex;
uniform sampler2D search_tex;
in vec2 texcoord;
in vec2 pixcoord;
in vec4 offset[3];
in vec4 dummy2;

out vec4 frag_color;

void main()
{
  frag_color = SMAABlendingWeightCalculationPS(texcoord, pixcoord, offset, edge_tex, area_tex, search_tex, ivec4(0));
}