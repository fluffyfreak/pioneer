// Copyright © 2008-2012 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "Player.h"
#include "Frame.h"
#include "Game.h"
#include "KeyBindings.h"
#include "Planet.h"
#include "Lang.h"
#include "Pi.h"
#include "SectorView.h"
#include "Serializer.h"
#include "ShipCpanel.h"
#include "Sound.h"
#include "SpaceStation.h"
#include "SpaceStationView.h"
#include "WorldView.h"

//Some player specific sounds
static Sound::Event s_soundUndercarriage;
static Sound::Event s_soundHyperdrive;

Player::Player(ShipType::Id shipId): Ship(shipId)
{
	SetController(new PlayerShipController());
	m_killCount = 0;
	m_knownKillCount = 0;
	EVA = false;
	EVABody = new EVAModel;
	EVABody->SetModel("sphere");
}

void Player::Save(Serializer::Writer &wr, Space *space)
{
	Ship::Save(wr, space);
	MarketAgent::Save(wr);
	wr.Int32(m_killCount);
	wr.Int32(m_knownKillCount);
}

void Player::Load(Serializer::Reader &rd, Space *space)
{
	Pi::player = this;
	Ship::Load(rd, space);
	MarketAgent::Load(rd);
	m_killCount = rd.Int32();
	m_knownKillCount = rd.Int32();
	EVA = false;
	// TODO: save EVABody status
	EVABody = new EVAModel;
	EVABody->SetModel("sphere");
}

//XXX remove + move to lua
void Player::OnHaveKilled(Body *guyWeKilled)
{
	if (guyWeKilled->IsType(Object::SHIP)) {
		printf("Well done. you killed some poor fucker\n");
		m_killCount++;
	}
}

//XXX perhaps remove this, the sound is very annoying
bool Player::OnDamage(Object *attacker, float kgDamage)
{
	bool r = Ship::OnDamage(attacker, kgDamage);
	if (!IsDead() && (GetPercentHull() < 25.0f)) {
		Sound::BodyMakeNoise(this, "warning", .5f);
	}
	return r;
}

//XXX handle killcounts in lua
void Player::SetDockedWith(SpaceStation *s, int port)
{
	Ship::SetDockedWith(s, port);
	if (s) {
		if (Pi::CombatRating(m_killCount) > Pi::CombatRating(m_knownKillCount)) {
			Pi::cpan->MsgLog()->ImportantMessage(Lang::PIONEERING_PILOTS_GUILD, Lang::RIGHT_ON_COMMANDER);
		}
		m_knownKillCount = m_killCount;

		Pi::SetView(Pi::spaceStationView);
	}
}

//XXX all ships should make this sound
bool Player::SetWheelState(bool down)
{
	bool did = Ship::SetWheelState(down);
	if (did) {
		s_soundUndercarriage.Play(down ? "UC_out" : "UC_in", 1.0f, 1.0f, 0);
	}
	return did;
}

bool Player::ToggleEVA()
{
	if (EVA){
		EVABody->DeActivate();
		EVA = false;

		return true;
	}
	else {
		if (GetFlightState() == Ship::LANDED)
		{
			matrix4x4d m;
			GetRotMatrix(m);
			vector3d pos = GetPosition() + m*vector3d(5, 0, 0);
			EVABody->Activate(GetFrame(), pos, m);
			EVA = true;
			return true;
		}
		else
			return false;
	}
}

//XXX all ships should make this sound
bool Player::FireMissile(int idx, Ship *target)
{
	if (!Ship::FireMissile(idx, target))
		return false;

	Sound::PlaySfx("Missile_launch", 1.0f, 1.0f, 0);
	return true;
}

//XXX do in lua, or use the alert concept for all ships
void Player::SetAlertState(Ship::AlertState as)
{
	Ship::AlertState prev = GetAlertState();

	switch (as) {
		case ALERT_NONE:
			if (prev != ALERT_NONE)
				Pi::cpan->MsgLog()->Message("", Lang::ALERT_CANCELLED);
			break;

		case ALERT_SHIP_NEARBY:
			if (prev == ALERT_NONE)
				Pi::cpan->MsgLog()->ImportantMessage("", Lang::SHIP_DETECTED_NEARBY);
			else
				Pi::cpan->MsgLog()->ImportantMessage("", Lang::DOWNGRADING_ALERT_STATUS);
			Sound::PlaySfx("OK");
			break;

		case ALERT_SHIP_FIRING:
			Pi::cpan->MsgLog()->ImportantMessage("", Lang::LASER_FIRE_DETECTED);
			Sound::PlaySfx("warning", 0.2f, 0.2f, 0);
			break;
	}

	Pi::cpan->SetAlertState(as);

	Ship::SetAlertState(as);
}

void Player::NotifyRemoved(const Body* const removedBody)
{
	if (GetNavTarget() == removedBody)
		SetNavTarget(0);

	else if (GetCombatTarget() == removedBody) {
		SetCombatTarget(0);

		if (!GetNavTarget() && removedBody->IsType(Object::SHIP))
			SetNavTarget(static_cast<const Ship*>(removedBody)->GetHyperspaceCloud());
	}

	Ship::NotifyRemoved(removedBody);
}

/* MarketAgent shite */
//XXX move to Player character .cpp
void Player::Bought(Equip::Type t)
{
	m_equipment.Add(t);
	UpdateEquipStats();
}

void Player::Sold(Equip::Type t)
{
	m_equipment.Remove(t, 1);
	UpdateEquipStats();
}

bool Player::CanBuy(Equip::Type t, bool verbose) const
{
	Equip::Slot slot = Equip::types[int(t)].slot;
	bool freespace = (m_equipment.FreeSpace(slot)!=0);
	bool freecapacity = (GetStats().free_capacity >= Equip::types[int(t)].mass);
	if (verbose) {
		if (!freespace) {
			Pi::Message(Lang::NO_FREE_SPACE_FOR_ITEM);
		}
		else if (!freecapacity) {
			Pi::Message(Lang::SHIP_IS_FULLY_LADEN);
		}
	}
	return (freespace && freecapacity);
}

bool Player::CanSell(Equip::Type t, bool verbose) const
{
	Equip::Slot slot = Equip::types[int(t)].slot;
	bool cansell = (m_equipment.Count(slot, t) > 0);
	if (verbose) {
		if (!cansell) {
			Pi::Message(stringf(Lang::YOU_DO_NOT_HAVE_ANY_X, formatarg("item", Equip::types[int(t)].name)));
		}
	}
	return cansell;
}

Sint64 Player::GetPrice(Equip::Type t) const
{
	if (Ship::GetDockedWith()) {
		return Ship::GetDockedWith()->GetPrice(t);
	} else {
		assert(0);
		return 0;
	}
}

EVAModel::EVAModel() : Ship(), eyeHeight(5), active(false) {
	SetMass(1000); // weighs a fricking tonne
}

void EVAModel::Activate(Frame* frame, vector3d pos, const matrix4x4d& rot) {
	active = true;
	SetFrame(frame);
	SetRotMatrix(rot);
	Body *astro = GetFrame()->m_astroBody;
	if (astro && astro->IsType(Object::TERRAINBODY)) {
		const TerrainBody *tb = static_cast<TerrainBody*>(astro);
		const vector3d npos = pos.Normalized();
		pos = npos * ( tb->GetTerrainHeight(npos) + eyeHeight);
	}
	SetPosition(pos);
	SetVelocity(vector3d(0,0,0));
	SetAngVelocity(vector3d(0,0,0));
	//PutOnSurfaceAt(pos);
	Enable();
	Pi::game->GetSpace()->AddBody(this);
}
void EVAModel::DeActivate() {
	active = false;
	Disable();
	Pi::game->GetSpace()->RemoveBody(this);
}

void EVAModel::StaticUpdate(const float timeStep) {
	if (!active)
		return;

	wantVel = wantAngVel = vector3d(0,0,0);

	if (KeyBindings::thrustLeft.IsActive()) wantVel.x += 1.0;
	if (KeyBindings::thrustRight.IsActive()) wantVel.x -= 1.0;
	if (KeyBindings::thrustUp.IsActive()) wantVel.y += 1.0;
	if (KeyBindings::thrustDown.IsActive()) wantVel.y -= 1.0;
	if (KeyBindings::thrustForward.IsActive()) wantVel.z += 1.0;
	if (KeyBindings::thrustBackwards.IsActive()) wantVel.z -= 1.0;

	if (KeyBindings::yawLeft.IsActive()) wantAngVel.y += 1.0;
	if (KeyBindings::yawRight.IsActive()) wantAngVel.y -= 1.0;
	if (KeyBindings::pitchUp.IsActive()) wantAngVel.x += 1.0;
	if (KeyBindings::pitchDown.IsActive()) wantAngVel.x -= 1.0;
	if (KeyBindings::rollLeft.IsActive()) wantAngVel.z += 1.0;
	if (KeyBindings::rollRight.IsActive()) wantAngVel.z -= 1.0;

	//if (KeyBindings::fastRotate.IsActive())
	//	wantAngVel *= 3.0;

	DynamicBody::StaticUpdate(timeStep);
}

void EVAModel::TimeStepUpdate(const float timeStep) {
	if (!active)
		return;
	// wandering the universe in a giant hamster-ball! Or is it just an
	// ill-modeled astronaut? Hard to tell.
	const double velDamping = 1;
	const double vertDamping = 3;
	const double angDamping = 1;
	const double straightening = 1;
	const double maxWalkVel = 50;
	const double maxTurnAngVel = 20;
	const double thrusterStrength = 5;
	const double springConstant = 200;

	const double f = GetMass();
	const double thrusterForce = f*thrusterStrength;
	const vector3d vel = GetVelocity();
	const vector3d angvel = GetAngVelocity();
	const vector3d pos = GetPosition();
	Body *astro = GetFrame()->m_astroBody;

	vector3d up, ahead, left, norm;

	matrix4x4d rot;
	GetRotMatrix(rot);
	const vector3d viewdir = rot*vector3d(0,0,-1);

	bool grounded = false;

	if (timeStep < 0.1) {
		if (astro && astro->IsType(Object::TERRAINBODY)) {
			const TerrainBody *tb = static_cast<TerrainBody*>(astro);
			const vector3d npos = pos.Normalized();
			const double h0 = tb->GetTerrainHeight(npos);
			double height = pos.Length() - h0;

			if (height < 1.1 * eyeHeight) {
				grounded = true;
				up = npos;
				if (up.Dot(viewdir) < 0.95)
					left = up.Cross(viewdir).Normalized();
				else
					left = up.Cross(rot*vector3d(-1,-1,-1));
				ahead = -up.Cross(left).Normalized();

				// estimate local surface normal
				const double dh1 = h0 - tb->GetTerrainHeight((pos+left).Normalized());
				const double dh2 = h0 - tb->GetTerrainHeight((pos+ahead).Normalized());
				norm = (up + left*dh1 + ahead*dh2).Normalized();
				if (norm.Dot(viewdir) < 0.95)
					left = norm.Cross(viewdir).Normalized();
				else
					left = norm.Cross(rot*vector3d(-1,-1,-1));
				ahead = -norm.Cross(left).Normalized();

				AddForce(f*norm*(springConstant*(eyeHeight-height)));
			}
		}
		if (!grounded) {
			ahead = viewdir;
			left = rot*vector3d(-1,0,0);
			up = rot*vector3d(0,1,0);
			norm = up;
		}

		const bool walking = grounded && (vel.Length() < maxWalkVel);
		AddForce(
				left * wantVel.x * ( walking ? f*15 : thrusterForce ) +
				norm * wantVel.y * ( walking ?
					( (wantVel.y > 0) ? f*50 : f*springConstant/2)
					: thrusterForce ) +
				ahead * wantVel.z * (walking ? f*30 : thrusterForce ));

		if (angvel.Length() < maxTurnAngVel)
			AddTorque(f*(
						left * wantAngVel.x +
						up * wantAngVel.y +
						viewdir * wantAngVel.z));

		if (grounded) {
			AddTorque(straightening*f*left*norm.Dot(viewdir));
			AddTorque(straightening*f*ahead*norm.Dot(rot*vector3d(1,0,0)));
			AddForce(vertDamping*-f*norm*norm.Dot(vel));
			AddTorque(angDamping*-f*angvel);
			if (wantVel.Length() <= 0.01)
				AddForce(-velDamping*f*vel);
		}
	} else {
		// time accelerated: physics breaks down, so we ignore it!
		if (astro && astro->IsType(Object::TERRAINBODY)) {
			const TerrainBody *tb = static_cast<TerrainBody*>(astro);
			const vector3d npos = pos.Normalized();
			const double h0 = tb->GetTerrainHeight(npos);
			double height = pos.Length() - h0;
			if (height < 2 * eyeHeight)
				SetPosition(npos*(h0+eyeHeight));
			SetVelocity(vector3d(0,0,0));
			SetAngVelocity(vector3d(0,0,0));
		}
	}

	DynamicBody::TimeStepUpdate(timeStep);
}

//void EVAModel::Render(const vector3d &viewCoords, const matrix4x4d &viewTransform) {
//	ModelBody::RenderLmrModel(viewCoords, viewTransform);
//}

//XXX ui stuff
void Player::OnEnterHyperspace()
{
	s_soundHyperdrive.Play("Hyperdrive_Jump");
	SetNavTarget(0);
	SetCombatTarget(0);

	Pi::worldView->HideTargetActions(); // hide the comms menu
	m_controller->SetFlightControlState(CONTROL_MANUAL); //could set CONTROL_HYPERDRIVE
	ClearThrusterState();
	Pi::game->WantHyperspace();
}

void Player::OnEnterSystem()
{
	m_controller->SetFlightControlState(CONTROL_MANUAL);
	//XXX don't call sectorview from here, use signals instead
	Pi::sectorView->ResetHyperspaceTarget();
}

//temporary targeting stuff
PlayerShipController *Player::GetPlayerController() const
{
	return static_cast<PlayerShipController*>(GetController());
}

Body *Player::GetCombatTarget() const
{
	return static_cast<PlayerShipController*>(m_controller)->GetCombatTarget();
}

Body *Player::GetNavTarget() const
{
	return static_cast<PlayerShipController*>(m_controller)->GetNavTarget();
}

Body *Player::GetSetSpeedTarget() const
{
	return static_cast<PlayerShipController*>(m_controller)->GetSetSpeedTarget();
}

void Player::SetCombatTarget(Body* const target, bool setSpeedTo)
{
	static_cast<PlayerShipController*>(m_controller)->SetCombatTarget(target, setSpeedTo);
	Pi::onPlayerChangeTarget.emit();
}

void Player::SetNavTarget(Body* const target, bool setSpeedTo)
{
	static_cast<PlayerShipController*>(m_controller)->SetNavTarget(target, setSpeedTo);
	Pi::onPlayerChangeTarget.emit();
}
//temporary targeting stuff ends

Ship::HyperjumpStatus Player::StartHyperspaceCountdown(const SystemPath &dest)
{
	HyperjumpStatus status = Ship::StartHyperspaceCountdown(dest);

	if (status == HYPERJUMP_OK)
		s_soundHyperdrive.Play("Hyperdrive_Charge");

	return status;
}

void Player::ResetHyperspaceCountdown()
{
	s_soundHyperdrive.Play("Hyperdrive_Abort");
	Ship::ResetHyperspaceCountdown();
}
