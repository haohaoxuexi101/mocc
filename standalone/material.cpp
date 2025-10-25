#include "material.hpp"

#include <stdexcept>

namespace standalone {

Material::Material(double total, double scatter, double nu_sig_f, double chi_val, double source)
    : sigma_t(total), sigma_s(scatter), nu_sigma_f(nu_sig_f), chi(chi_val), external_source(source) {}

double Material::isotropic_source(double scalar_flux) const noexcept {
    return external_source + sigma_s * scalar_flux;
}

std::size_t MaterialLibrary::add_material(const std::string &name, const Material &material) {
    auto iter = lookup_.find(name);
    if (iter != lookup_.end()) {
        throw std::runtime_error("material already exists: " + name);
    }
    std::size_t idx = materials_.size();
    materials_.push_back(material);
    lookup_.emplace(name, idx);
    return idx;
}

std::size_t MaterialLibrary::id(const std::string &name) const {
    auto iter = lookup_.find(name);
    if (iter == lookup_.end()) {
        throw std::runtime_error("unknown material: " + name);
    }
    return iter->second;
}

const Material &MaterialLibrary::material(std::size_t idx) const {
    if (idx >= materials_.size()) {
        throw std::out_of_range("material index out of range");
    }
    return materials_[idx];
}

const Material &MaterialLibrary::material(const std::string &name) const {
    return material(id(name));
}

} // namespace standalone
