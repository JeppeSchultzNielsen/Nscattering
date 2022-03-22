#ifndef PARTICLETYPE_H
#define PARTICLETYPE_H

#include <ausa/eloss/Ion.h>

class ParticleType {
public:
    explicit ParticleType(const std::string &name);

    static ParticleType constructUnknown() {
        auto ret = ParticleType{};
        ret.unknown = true;
        return ret;
    }

    double mass{};
    std::string name{};
    int A{}, Z{};
    bool unknown{};

    ParticleType() = default;
};

ParticleType::ParticleType(const std::string &name) {
    AUSA::EnergyLoss::Ion ion{name};
    mass = ion.getMass();
    this->name = ion.getName();
    Z = ion.getZ();
    A = ion.getA();
    unknown = false;
}

#endif //PARTICLETYPE_H
