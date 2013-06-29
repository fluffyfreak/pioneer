// Copyright © 2008-2013 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef _SHIPTYPE_H
#define _SHIPTYPE_H

#include "libs.h"
#include "vector3.h"
#include "EquipType.h"
#include <vector>
#include <map>

struct GunMountData {	// fixed mount, front or rear
	vector3d pos;
	vector3d dir;		// not necessarily directly front or rear...
	Sint32 size;		// maximum size of weapon to mount  
	std::string name;
};

struct TurretData : public GunMountData {
	double extent;		// maximum angle from dir in radians
	double accel;
	double maxspeed;	// max training speed & accel in radians
};

struct ShipType {
	enum Thruster { // <enum scope='ShipType' name=ShipTypeThruster prefix=THRUSTER_ public>
		THRUSTER_REVERSE,
		THRUSTER_FORWARD,
		THRUSTER_UP,
		THRUSTER_DOWN,
		THRUSTER_LEFT,
		THRUSTER_RIGHT,
		THRUSTER_MAX // <enum skip>
	};
	enum DualLaserOrientation { // <enum scope='ShipType' name='DualLaserOrientation' prefix='DUAL_LASERS_' public>
		DUAL_LASERS_HORIZONTAL,
		DUAL_LASERS_VERTICAL
	};
	enum Tag { // <enum scope='ShipType' name=ShipTypeTag prefix=TAG_ public>
		TAG_NONE,
		TAG_SHIP,
		TAG_STATIC_SHIP,
		TAG_MISSILE,
		TAG_MAX // <enum skip>
	};
	typedef std::string Id;

	////////
	Tag tag;
	Id id;
	std::string name;
	std::string modelName;
	float linThrust[THRUSTER_MAX];
	float angThrust;
	std::vector<GunMountData> gunMount;
	std::vector<TurretData> turret;
	int equipSlotCapacity[Equip::SLOT_MAX];
	int capacity; // tonnes
	int hullMass;
	float effectiveExhaustVelocity; // velocity at which the propellant escapes the engines
	int fuelTankMass; //full fuel tank mass, on top of hullMass
	int baseprice;
	Equip::Type hyperdrive;
	vector3d cameraOffset;
	int minCrew, maxCrew; // XXX really only for Lua, but needs to be declared in the ship def
	///////

	static std::string POLICE;
	static std::string MISSILE_GUIDED;
	static std::string MISSILE_NAVAL;
	static std::string MISSILE_SMART;
	static std::string MISSILE_UNGUIDED;

	static std::map<Id, ShipType> types;
	static std::vector<Id> player_ships;
	static std::vector<Id> static_ships;
	static std::vector<Id> missile_ships;

	static std::vector<Id> playable_atmospheric_ships;

	static void Init();
	static const ShipType *Get(const char *name) {
		std::map<Id, ShipType>::iterator t = types.find(name);
		if (t == types.end()) return 0;
		else return &(*t).second;
	}
};

#endif /* _SHIPTYPE_H */
