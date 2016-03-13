// Copyright © 2008-2016 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "Terrain.h"
#include "perlin.h"
#include "Pi.h"
#include "FileSystem.h"
#include "FloatComparison.h"

// static instancer. selects the best height and color classes for the body
Terrain *Terrain::InstanceTerrain(const SystemBody *body)
{
	static const std::string s_dummy;
	std::string jsonHeightFilename = s_dummy;
	std::string jsonColourFilename = s_dummy;

	Random rand(body->GetSeed());
	GeneratorInstancer gi = InstanceGenerator<TerrainHeightJSON,TerrainColourJSON>;

	// special case for heightmaps
	// XXX this is terrible but will do for now until we get a unified
	// heightmap setup. if you add another height fractal, remember to change
	// the check in CustomSystem::l_height_map
	if (!body->GetJSONFilename().empty()) 
	{
		jsonHeightFilename = body->GetJSONFilename();
		jsonColourFilename = body->GetJSONFilename();
	}
	else 
	{
		switch (body->GetType()) 
		{
			case SystemBody::TYPE_BROWN_DWARF:
				jsonHeightFilename = "ellipsoid.json";
				jsonColourFilename = "StarBrownDwarf.json";
				break;

			case SystemBody::TYPE_WHITE_DWARF:
				jsonHeightFilename = "ellipsoid.json";
				jsonColourFilename = "StarWhiteDwarf.json";
				break;

			case SystemBody::TYPE_STAR_M:
			case SystemBody::TYPE_STAR_M_GIANT:
			case SystemBody::TYPE_STAR_M_SUPER_GIANT:
			case SystemBody::TYPE_STAR_M_HYPER_GIANT: {
				static const std::pair<std::string, std::string> choices[] = {
					std::make_pair<std::string, std::string>("ellipsoid.json","StarK.json"),
					std::make_pair<std::string, std::string>("ellipsoid.json","StarK.json"),
					std::make_pair<std::string, std::string>("ellipsoid.json","StarK.json"),
					std::make_pair<std::string, std::string>("ellipsoid.json","StarG.json")
				};
				const Uint32 choice = rand.Int32(COUNTOF(choices)); jsonHeightFilename = choices[choice].first; jsonColourFilename = choices[choice].second;
				break;
			}

			case SystemBody::TYPE_STAR_K:
			case SystemBody::TYPE_STAR_K_GIANT:
			case SystemBody::TYPE_STAR_K_SUPER_GIANT:
			case SystemBody::TYPE_STAR_K_HYPER_GIANT: {
				static const std::pair<std::string, std::string> choices[] = {
					std::make_pair<std::string, std::string>("ellipsoid.json","StarK.json"),
					std::make_pair<std::string, std::string>("ellipsoid.json","StarK.json"),
					std::make_pair<std::string, std::string>("ellipsoid.json","StarK.json"),
					std::make_pair<std::string, std::string>("ellipsoid.json","StarG.json")
				};
				const Uint32 choice = rand.Int32(COUNTOF(choices)); jsonHeightFilename = choices[choice].first; jsonColourFilename = choices[choice].second;
				break;
			}

			case SystemBody::TYPE_STAR_G:
			case SystemBody::TYPE_STAR_G_GIANT:
			case SystemBody::TYPE_STAR_G_SUPER_GIANT:
			case SystemBody::TYPE_STAR_G_HYPER_GIANT: {
				static const std::pair<std::string, std::string> choices[] = {
					std::make_pair<std::string, std::string>("ellipsoid.json","StarWhiteDwarf.json"),
					std::make_pair<std::string, std::string>("ellipsoid.json","StarG.json")
				};
				const Uint32 choice = rand.Int32(COUNTOF(choices)); jsonHeightFilename = choices[choice].first; jsonColourFilename = choices[choice].second;
				break;
			}

			case SystemBody::TYPE_STAR_F:
			case SystemBody::TYPE_STAR_F_GIANT:
			case SystemBody::TYPE_STAR_F_HYPER_GIANT:
			case SystemBody::TYPE_STAR_F_SUPER_GIANT:
			case SystemBody::TYPE_STAR_A:
			case SystemBody::TYPE_STAR_A_GIANT:
			case SystemBody::TYPE_STAR_A_HYPER_GIANT:
			case SystemBody::TYPE_STAR_A_SUPER_GIANT:
			case SystemBody::TYPE_STAR_B:
			case SystemBody::TYPE_STAR_B_GIANT:
			case SystemBody::TYPE_STAR_B_SUPER_GIANT:
			case SystemBody::TYPE_STAR_B_WF:
			case SystemBody::TYPE_STAR_O:
			case SystemBody::TYPE_STAR_O_GIANT:
			case SystemBody::TYPE_STAR_O_HYPER_GIANT:
			case SystemBody::TYPE_STAR_O_SUPER_GIANT:
			case SystemBody::TYPE_STAR_O_WF:
				jsonHeightFilename = "ellipsoid.json";
				jsonColourFilename = "Solid.json";
			break;

			case SystemBody::TYPE_PLANET_GAS_GIANT: {
				static const std::pair<std::string, std::string> choices[] = {
					std::make_pair<std::string, std::string>("flat.json","GGJupiter.json"),
					std::make_pair<std::string, std::string>("flat.json","GGSaturn.json"),
					std::make_pair<std::string, std::string>("flat.json","GGSaturn2.json"),
					std::make_pair<std::string, std::string>("flat.json","GGNeptune.json"),
					std::make_pair<std::string, std::string>("flat.json","GGNeptune2.json"),
					std::make_pair<std::string, std::string>("flat.json","GGUranus.json"),
					std::make_pair<std::string, std::string>("flat.json","GGSaturn.json")
				};
				const Uint32 choice = rand.Int32(COUNTOF(choices)); jsonHeightFilename = choices[choice].first; jsonColourFilename = choices[choice].second;
				break;
			}

			case SystemBody::TYPE_PLANET_ASTEROID: {
				static const std::pair<std::string, std::string> choices[] = {
					std::make_pair<std::string, std::string>("asteroid.json","Asteroid.json"),
					std::make_pair<std::string, std::string>("asteroid2.json","Asteroid.json"),
					std::make_pair<std::string, std::string>("asteroid3.json","Asteroid.json"),
					std::make_pair<std::string, std::string>("asteroid4.json","Asteroid.json"),
					std::make_pair<std::string, std::string>("asteroid.json","Rock.json"),
					std::make_pair<std::string, std::string>("asteroid2.json","BandedRock.json"),
					std::make_pair<std::string, std::string>("asteroid3.json","Rock.json"),
					std::make_pair<std::string, std::string>("asteroid4.json","BandedRock.json")
				};
				const Uint32 choice = rand.Int32(COUNTOF(choices)); jsonHeightFilename = choices[choice].first; jsonColourFilename = choices[choice].second;
				break;
			}

			case SystemBody::TYPE_PLANET_TERRESTRIAL: {

				//Over-ride:
				//jsonHeightFilename = "asteroid3.json"; jsonColourFilename = "Rock.json";
				//break;
				// Earth-like world

				if ((body->GetLifeAsFixed() > fixed(7,10)) && (body->GetVolatileGasAsFixed() > fixed(2,10))) {
					// There would be no life on the surface without atmosphere

					if (body->GetAverageTemp() > 240) {
						static const std::pair<std::string, std::string> choices[] = {
							std::make_pair<std::string, std::string>("hillsridged.json","EarthLike.json"),
							std::make_pair<std::string, std::string>("hillsrivers.json","EarthLike.json"),
							std::make_pair<std::string, std::string>("hillsdunes.json","EarthLike.json"),
							std::make_pair<std::string, std::string>("mountainsridged.json","EarthLike.json"),
							std::make_pair<std::string, std::string>("mountainsnormal.json","EarthLike.json"),
							std::make_pair<std::string, std::string>("mountainsrivers.json","EarthLike.json"),
							std::make_pair<std::string, std::string>("mountainsvolcano.json","EarthLike.json"),
							std::make_pair<std::string, std::string>("mountainsriversvolcano.json","EarthLike.json")
						};
						const Uint32 choice = rand.Int32(COUNTOF(choices)); jsonHeightFilename = choices[choice].first; jsonColourFilename = choices[choice].second;
						break;
					}

					// desert-ice planets
					static const std::pair<std::string, std::string> choices[] = {
						std::make_pair<std::string, std::string>("hillsridged.json","Desert.json"),
						std::make_pair<std::string, std::string>("hillsrivers.json","Desert.json"),
						std::make_pair<std::string, std::string>("hillsdunes.json","Desert.json"),
						std::make_pair<std::string, std::string>("mountainsridged.json","Desert.json"),
						std::make_pair<std::string, std::string>("mountainsnormal.json","Desert.json"),
						std::make_pair<std::string, std::string>("mountainsrivers.json","Desert.json"),
						std::make_pair<std::string, std::string>("mountainsvolcano.json","Desert.json"),
						std::make_pair<std::string, std::string>("mountainsriversvolcano.json","Desert.json"),
						std::make_pair<std::string, std::string>("barrenrock.json","Desert.json"),
						std::make_pair<std::string, std::string>("barrenrock2.json","Desert.json")
					};
					const Uint32 choice = rand.Int32(COUNTOF(choices)); jsonHeightFilename = choices[choice].first; jsonColourFilename = choices[choice].second;
					break;
				}

				// Harsh, habitable world
				if ((body->GetVolatileGasAsFixed() > fixed(2,10)) && (body->GetLifeAsFixed() > fixed(4,10)) ) {

					if (body->GetAverageTemp() > 240) {
						static const std::pair<std::string, std::string> choices[] = {
							std::make_pair<std::string, std::string>("hillsridged.json","TFGood.json"),
							std::make_pair<std::string, std::string>("hillsrivers.json","TFGood.json"),
							std::make_pair<std::string, std::string>("hillsdunes.json","TFGood.json"),
							std::make_pair<std::string, std::string>("hillsnormal.json","TFGood.json"),
							std::make_pair<std::string, std::string>("mountainsnormal.json","TFGood.json"),
							std::make_pair<std::string, std::string>("mountainsridged.json","TFGood.json"),
							std::make_pair<std::string, std::string>("mountainsvolcano.json","TFGood.json"),
							std::make_pair<std::string, std::string>("mountainsriversvolcano.json","TFGood.json"),
							std::make_pair<std::string, std::string>("mountainsrivers.json","TFGood.json"),
							std::make_pair<std::string, std::string>("ruggeddesert.json","TFGood.json"),
							std::make_pair<std::string, std::string>("barrenrock.json","TFGood.json"),
							std::make_pair<std::string, std::string>("barrenrock2.json","TFGood.json")
						};
						const Uint32 choice = rand.Int32(COUNTOF(choices)); jsonHeightFilename = choices[choice].first; jsonColourFilename = choices[choice].second;
						break;
					}

					// ice planets
					static const std::pair<std::string, std::string> choices[] = {
						std::make_pair<std::string, std::string>("hillsridged.json","Ice.json"),
						std::make_pair<std::string, std::string>("hillsrivers.json","Ice.json"),
						std::make_pair<std::string, std::string>("hillsdunes.json","Ice.json"),
						std::make_pair<std::string, std::string>("hillsnormal.json","Ice.json"),
						std::make_pair<std::string, std::string>("mountainsnormal.json","Ice.json"),
						std::make_pair<std::string, std::string>("mountainsridged.json","Ice.json"),
						std::make_pair<std::string, std::string>("mountainsvolcano.json","Ice.json"),
						std::make_pair<std::string, std::string>("mountainsriversvolcano.json","Ice.json"),
						std::make_pair<std::string, std::string>("mountainsrivers.json","Ice.json"),
						std::make_pair<std::string, std::string>("ruggeddesert.json","Ice.json"),
						std::make_pair<std::string, std::string>("barrenrock.json","Ice.json"),
						std::make_pair<std::string, std::string>("barrenrock2.json","Ice.json"),
						std::make_pair<std::string, std::string>("barrenrock3.json","Ice.json")
					};
					const Uint32 choice = rand.Int32(COUNTOF(choices)); jsonHeightFilename = choices[choice].first; jsonColourFilename = choices[choice].second;
					break;
				}

				// Marginally habitable world/ verging on mars like :)
				else if ((body->GetVolatileGasAsFixed() > fixed(1,10)) && (body->GetLifeAsFixed() > fixed(1,10)) ) {

					if (body->GetAverageTemp() > 240) {
						static const std::pair<std::string, std::string> choices[] = {
							std::make_pair<std::string, std::string>("hillsridged.json","TFPoor.json"),
							std::make_pair<std::string, std::string>("hillsrivers.json","TFPoor.json"),
							std::make_pair<std::string, std::string>("hillsdunes.json","TFPoor.json"),
							std::make_pair<std::string, std::string>("hillsnormal.json","TFPoor.json"),
							std::make_pair<std::string, std::string>("mountainsnormal.json","TFPoor.json"),
							std::make_pair<std::string, std::string>("mountainsridged.json","TFPoor.json"),
							std::make_pair<std::string, std::string>("mountainsvolcano.json","TFPoor.json"),
							std::make_pair<std::string, std::string>("mountainsriversvolcano.json","TFPoor.json"),
							std::make_pair<std::string, std::string>("mountainsrivers.json","TFPoor.json"),
							std::make_pair<std::string, std::string>("ruggeddesert.json","TFPoor.json"),
							std::make_pair<std::string, std::string>("barrenrock.json","TFPoor.json"),
							std::make_pair<std::string, std::string>("barrenrock2.json","TFPoor.json"),
							std::make_pair<std::string, std::string>("barrenrock3.json","TFPoor.json")
						};
						const Uint32 choice = rand.Int32(COUNTOF(choices)); jsonHeightFilename = choices[choice].first; jsonColourFilename = choices[choice].second;
						break;
					}

					// ice planets
					static const std::pair<std::string, std::string> choices[] = {
						std::make_pair<std::string, std::string>("hillsridged.json","Ice.json"),
						std::make_pair<std::string, std::string>("hillsrivers.json","Ice.json"),
						std::make_pair<std::string, std::string>("hillsdunes.json","Ice.json"),
						std::make_pair<std::string, std::string>("hillsnormal.json","Ice.json"),
						std::make_pair<std::string, std::string>("mountainsnormal.json","Ice.json"),
						std::make_pair<std::string, std::string>("mountainsridged.json","Ice.json"),
						std::make_pair<std::string, std::string>("mountainsvolcano.json","Ice.json"),
						std::make_pair<std::string, std::string>("mountainsriversvolcano.json","Ice.json"),
						std::make_pair<std::string, std::string>("mountainsrivers.json","Ice.json"),
						std::make_pair<std::string, std::string>("ruggeddesert.json","Ice.json"),
						std::make_pair<std::string, std::string>("barrenrock.json","Ice.json"),
						std::make_pair<std::string, std::string>("barrenrock2.json","Ice.json"),
						std::make_pair<std::string, std::string>("barrenrock3.json","Ice.json")
					};
					const Uint32 choice = rand.Int32(COUNTOF(choices)); jsonHeightFilename = choices[choice].first; jsonColourFilename = choices[choice].second;
					break;
				}

				// Desert-like world, Mars -like.
				if ((body->GetVolatileLiquidAsFixed() < fixed(1,10)) && (body->GetVolatileGasAsFixed() > fixed(1,5))) {
					static const std::pair<std::string, std::string> choices[] = {
						std::make_pair<std::string, std::string>("hillsdunes.json","Desert.json"),
						std::make_pair<std::string, std::string>("watersolid.json","Desert.json"),
						std::make_pair<std::string, std::string>("ruggeddesert.json","Desert.json"),
						std::make_pair<std::string, std::string>("ruggedlava.json","Desert.json"),
						std::make_pair<std::string, std::string>("mountainsvolcano.json","Desert.json"),
						std::make_pair<std::string, std::string>("mountainsriversvolcano.json","Desert.json"),
						std::make_pair<std::string, std::string>("barrenrock.json","Desert.json"),
						std::make_pair<std::string, std::string>("barrenrock2.json","Desert.json")
					};
					const Uint32 choice = rand.Int32(COUNTOF(choices)); jsonHeightFilename = choices[choice].first; jsonColourFilename = choices[choice].second;
					break;
				}

				// Frozen world
				if ((body->GetVolatileIcesAsFixed() > fixed(8,10)) &&  (body->GetAverageTemp() < 250)) {
					static const std::pair<std::string, std::string> choices[] = {
						std::make_pair<std::string, std::string>("hillsdunes.json","Ice.json"),
						std::make_pair<std::string, std::string>("hillscraters.json","Ice.json"),
						std::make_pair<std::string, std::string>("mountainscraters.json","Ice.json"),
						std::make_pair<std::string, std::string>("watersolid.json","Ice.json"),
						std::make_pair<std::string, std::string>("watersolidcanyons.json","Ice.json"),
						std::make_pair<std::string, std::string>("ruggeddesert.json","Ice.json"),
						std::make_pair<std::string, std::string>("barrenrock.json","Ice.json"),
						std::make_pair<std::string, std::string>("barrenrock2.json","Ice.json"),
						std::make_pair<std::string, std::string>("barrenrock3.json","Ice.json")
					};
					const Uint32 choice = rand.Int32(COUNTOF(choices)); jsonHeightFilename = choices[choice].first; jsonColourFilename = choices[choice].second;
					break;
				}

				// Volcanic world
				if (body->GetVolcanicityAsFixed() > fixed(7,10)) {

					if (body->GetLifeAsFixed() > fixed(5,10)) {	// life on a volcanic world ;)
						jsonHeightFilename = "ruggedlava.json";
						jsonColourFilename = "TFGood.json";
					} else if (body->GetLifeAsFixed() > fixed(2,10)) {
						jsonHeightFilename = "ruggedlava.json";
						jsonColourFilename = "TFPoor.json";
					} else {
						jsonHeightFilename = "ruggedlava.json";
						jsonColourFilename = "Volcanic.json";
					}
					break;
				}

				//Below might not be needed.
				//Alien life world:
				if (body->GetLifeAsFixed() > fixed(1,10))  {
					static const std::pair<std::string, std::string> choices[] = {
						std::make_pair<std::string, std::string>("hillsdunes.json","TFPoor.json"),
						std::make_pair<std::string, std::string>("hillsridged.json","TFPoor.json"),
						std::make_pair<std::string, std::string>("hillsrivers.json","TFPoor.json"),
						std::make_pair<std::string, std::string>("mountainsnormal.json","TFPoor.json"),
						std::make_pair<std::string, std::string>("mountainsridged.json","TFPoor.json"),
						std::make_pair<std::string, std::string>("mountainsvolcano.json","TFPoor.json"),
						std::make_pair<std::string, std::string>("mountainsriversvolcano.json","TFPoor.json"),
						std::make_pair<std::string, std::string>("mountainsrivers.json","TFPoor.json"),
						std::make_pair<std::string, std::string>("watersolid.json","TFPoor.json"),
						std::make_pair<std::string, std::string>("ruggedlava.json","TFPoor.json"),
						std::make_pair<std::string, std::string>("ruggeddesert.json","TFPoor.json"),
						std::make_pair<std::string, std::string>("barrenrock.json","Ice.json"),
						std::make_pair<std::string, std::string>("barrenrock2.json","Ice.json"),
						std::make_pair<std::string, std::string>("barrenrock3.json","Ice.json")
					};
					const Uint32 choice = rand.Int32(COUNTOF(choices)); jsonHeightFilename = choices[choice].first; jsonColourFilename = choices[choice].second;
					break;
				};

				if (body->GetVolatileGasAsFixed() > fixed(1,10)) {
					static const std::pair<std::string, std::string> choices[] = {
						std::make_pair<std::string, std::string>("hillsnormal.json","Rock.json"),
						std::make_pair<std::string, std::string>("mountainsnormal.json","Rock.json"),
						std::make_pair<std::string, std::string>("ruggeddesert.json","Rock.json"),
						std::make_pair<std::string, std::string>("barrenrock.json","Rock.json"),
						std::make_pair<std::string, std::string>("barrenrock2.json","Rock.json"),
						std::make_pair<std::string, std::string>("barrenrock3.json","Rock.json")
					};
					const Uint32 choice = rand.Int32(COUNTOF(choices)); jsonHeightFilename = choices[choice].first; jsonColourFilename = choices[choice].second;
					break;
				}

				static const std::pair<std::string, std::string> choices[] = {
					std::make_pair<std::string, std::string>("hillscraters2.json","Rock.json"),
					std::make_pair<std::string, std::string>("mountainscraters2.json","Rock.json"),
					std::make_pair<std::string, std::string>("barrenrock3.json","Rock.json")
				};
				const Uint32 choice = rand.Int32(COUNTOF(choices)); jsonHeightFilename = choices[choice].first; jsonColourFilename = choices[choice].second;
				break;
			}

			default:
				jsonHeightFilename = "flat.json";
				jsonColourFilename = "Solid.json";
				break;
		}
	}

	return gi(body, jsonHeightFilename, jsonColourFilename);
}

Terrain::Terrain(const SystemBody *body) : m_seed(body->GetSeed()), m_rand(body->GetSeed()), m_minBody(body) 
{
	switch (Pi::detail.textures) {
		case 0: textures = false;
			m_fracnum = 2;break;
		default:
		case 1: textures = true;
			m_fracnum = 0;break;
	}

	switch (Pi::detail.fracmult) {
		case 0: m_fracmult = 100;break;
		case 1: m_fracmult = 10;break;
		case 2: m_fracmult = 1;break;
		case 3: m_fracmult = 0.5;break;
		default:
		case 4: m_fracmult = 0.1;break;
	}

	m_sealevel = Clamp(body->GetVolatileLiquid(), 0.0, 1.0);
	m_icyness  = Clamp(body->GetVolatileIces(), 0.0, 1.0);
	m_volcanic = Clamp(body->GetVolcanicity(), 0.0, 1.0); // height scales with volcanicity as well
	m_surfaceEffects = 0;

	const double rad = m_minBody.m_radius;

	// calculate max height
	m_maxHeightInMeters = std::max(100.0, (9000.0*rad*rad*(m_volcanic+0.5)) / (body->GetMass() * 6.64e-12));
	if (!is_finite(m_maxHeightInMeters)) 
		m_maxHeightInMeters = rad * 0.5;
	//             ^^^^ max mountain height for earth-like planet (same mass, radius)
	// and then in sphere normalized jizz
	m_maxHeight = std::min(1.0, m_maxHeightInMeters / rad);
	//Output("%s: max terrain height: %fm [%f]\n", m_minBody.name.c_str(), m_maxHeightInMeters, m_maxHeight);
	m_invMaxHeight = 1.0 / m_maxHeight;
	m_planetRadius = rad;
	m_planetEarthRadii = rad / EARTH_RADIUS;

	// Pick some colors, mainly reds and greens
	for (int i=0; i<int(COUNTOF(m_entropy)); i++) m_entropy[i] = m_rand.Double();
	for (int i=0; i<int(COUNTOF(m_rockColor)); i++) {
		double r,g,b;
		r = m_rand.Double(0.3, 1.0);
		g = m_rand.Double(0.3, r);
		b = m_rand.Double(0.3, g);
		r = std::max(b, r * body->GetMetallicity());
		g = std::max(b, g * body->GetMetallicity());
		m_rockColor[i] = vector3d(r, g, b);
	}

	// Pick some darker colours mainly reds and greens
	for (int i=0; i<int(COUNTOF(m_entropy)); i++) m_entropy[i] = m_rand.Double();
	for (int i=0; i<int(COUNTOF(m_darkrockColor)); i++) {
		double r,g,b;
		r = m_rand.Double(0.05, 0.3);
		g = m_rand.Double(0.05, r);
		b = m_rand.Double(0.05, g);
		r = std::max(b, r * body->GetMetallicity());
		g = std::max(b, g * body->GetMetallicity());
		m_darkrockColor[i] = vector3d(r, g, b);
	}

	// grey colours, in case you simply must have a grey colour on a world with high metallicity
	for (int i=0; i<int(COUNTOF(m_entropy)); i++) m_entropy[i] = m_rand.Double();
	for (int i=0; i<int(COUNTOF(m_greyrockColor)); i++) {
		double g;
		g = m_rand.Double(0.3, 0.9);
		m_greyrockColor[i] = vector3d(g, g, g);
	}

	// Pick some plant colours, mainly greens
	// TODO take star class into account
	for (int i=0; i<int(COUNTOF(m_entropy)); i++) m_entropy[i] = m_rand.Double();
	for (int i=0; i<int(COUNTOF(m_plantColor)); i++) {
		double r,g,b;
		g = m_rand.Double(0.3, 1.0);
		r = m_rand.Double(0.3, g);
		b = m_rand.Double(0.2, r);
		g = std::max(r, g * body->GetLife());
		b *= (1.0-body->GetLife());
		m_plantColor[i] = vector3d(r, g, b);
	}

	// Pick some darker plant colours mainly greens
	// TODO take star class into account
	for (int i=0; i<int(COUNTOF(m_entropy)); i++) m_entropy[i] = m_rand.Double();
	for (int i=0; i<int(COUNTOF(m_darkplantColor)); i++) {
		double r,g,b;
		g = m_rand.Double(0.05, 0.3);
		r = m_rand.Double(0.00, g);
		b = m_rand.Double(0.00, r);
		g = std::max(r, g * body->GetLife());
		b *= (1.0-body->GetLife());
		m_darkplantColor[i] = vector3d(r, g, b);
	}

	// Pick some sand colours, mainly yellow
	// TODO let some planetary value scale this colour
	for (int i=0; i<int(COUNTOF(m_entropy)); i++) m_entropy[i] = m_rand.Double();
	for (int i=0; i<int(COUNTOF(m_sandColor)); i++) {
		double r,g,b;
		r = m_rand.Double(0.6, 1.0);
		g = m_rand.Double(0.6, r);
		//b = m_rand.Double(0.0, g/2.0);
		b = 0;
		m_sandColor[i] = vector3d(r, g, b);
	}

	// Pick some darker sand colours mainly yellow
	// TODO let some planetary value scale this colour
	for (int i=0; i<int(COUNTOF(m_entropy)); i++) m_entropy[i] = m_rand.Double();
	for (int i=0; i<int(COUNTOF(m_darksandColor)); i++) {
		double r,g,b;
		r = m_rand.Double(0.05, 0.6);
		g = m_rand.Double(0.00, r);
		//b = m_rand.Double(0.00, g/2.0);
		b = 0;
		m_darksandColor[i] = vector3d(r, g, b);
	}

	// Pick some dirt colours, mainly red/brown
	// TODO let some planetary value scale this colour
	for (int i=0; i<int(COUNTOF(m_entropy)); i++) m_entropy[i] = m_rand.Double();
	for (int i=0; i<int(COUNTOF(m_dirtColor)); i++) {
		double r,g,b;
		r = m_rand.Double(0.3, 0.7);
		g = m_rand.Double(r-0.1, 0.75);
		b = m_rand.Double(0.0, r/2.0);
		m_dirtColor[i] = vector3d(r, g, b);
	}

	// Pick some darker dirt colours mainly red/brown
	// TODO let some planetary value scale this colour
	for (int i=0; i<int(COUNTOF(m_entropy)); i++) m_entropy[i] = m_rand.Double();
	for (int i=0; i<int(COUNTOF(m_darkdirtColor)); i++) {
		double r,g,b;
		r = m_rand.Double(0.05, 0.3);
		g = m_rand.Double(r-0.05, 0.35);
		b = m_rand.Double(0.0, r/2.0);
		m_darkdirtColor[i] = vector3d(r, g, b);
	}

	// These are used for gas giant colours, they are more m_random and *should* really use volatileGasses - TODO
	for (int i=0; i<int(COUNTOF(m_entropy)); i++) m_entropy[i] = m_rand.Double();
	for (int i=0; i<int(COUNTOF(m_gglightColor)); i++) {
		double r,g,b;
		r = m_rand.Double(0.0, 0.5);
		g = m_rand.Double(0.0, 0.5);
		b = m_rand.Double(0.0, 0.5);
		m_gglightColor[i] = vector3d(r, g, b);
	}
	//darker gas giant colours, more reds and greens
	for (int i=0; i<int(COUNTOF(m_entropy)); i++) m_entropy[i] = m_rand.Double();
	for (int i=0; i<int(COUNTOF(m_ggdarkColor)); i++) {
		double r,g,b;
		r = m_rand.Double(0.0, 0.3);
		g = m_rand.Double(0.0, r);
		b = m_rand.Double(0.0, std::min(r, g));
		m_ggdarkColor[i] = vector3d(r, g, b);
	}
}

Terrain::~Terrain()
{
}


/**
 * Feature width means roughly one perlin noise blob or grain.
 * This will end up being one hill, mountain or continent, roughly.
 */
void Terrain::SetFracDef(const unsigned int index, const double featureHeightMeters, const double featureWidthMeters, const double smallestOctaveMeters)
{
	assert(index>=0 && index<MAX_FRACDEFS);
	// feature
	m_fracdef[index].amplitude = featureHeightMeters / (m_maxHeight * m_planetRadius);
	m_fracdef[index].frequency = m_planetRadius / featureWidthMeters;
	m_fracdef[index].octaves = std::max(1, int(ceil(log(featureWidthMeters / smallestOctaveMeters) / log(2.0))));
	m_fracdef[index].lacunarity = 2.0;
	//Output("%d octaves\n", m_fracdef[index].octaves); //print
}

void Terrain::DebugDump() const
{
	Output("Terrain state dump:\n");
	Output("  Height fractal: %s\n", GetHeightFractalName());
	Output("  Color fractal: %s\n", GetColorFractalName());
	Output("  Detail: fracnum %d  fracmult %f  textures %s\n", m_fracnum, m_fracmult, textures ? "true" : "false");
	Output("  Config: DetailPlanets %d   FractalMultiple %d  Textures  %d\n", Pi::config->Int("DetailPlanets"), Pi::config->Int("FractalMultiple"), Pi::config->Int("Textures"));
	Output("  Seed: %d\n", m_seed);
	Output("  Body: %s [%d,%d,%d,%u,%u]\n", m_minBody.m_name.c_str(), m_minBody.m_path.sectorX, m_minBody.m_path.sectorY, m_minBody.m_path.sectorZ, m_minBody.m_path.systemIndex, m_minBody.m_path.bodyIndex);
	Output("  Aspect Ratio: %g\n", m_minBody.m_aspectRatio);
	Output("  Fracdefs:\n");
	for (int i = 0; i < 10; i++) {
		Output("    %d: amp %f  freq %f  lac %f  oct %d\n", i, m_fracdef[i].amplitude, m_fracdef[i].frequency, m_fracdef[i].lacunarity, m_fracdef[i].octaves);
	}
}
