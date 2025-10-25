#pragma once

#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

namespace standalone {

struct Material {
    double sigma_t{0.0};
    double sigma_s{0.0};
    double nu_sigma_f{0.0};
    double chi{1.0};
    double external_source{0.0};

    Material() = default;
    Material(double total, double scatter, double nu_sig_f, double chi_val, double source);

    double removal() const noexcept { return sigma_t - sigma_s; }
    double isotropic_source(double scalar_flux) const noexcept;
};

class MaterialLibrary {
public:
    std::size_t add_material(const std::string &name, const Material &material);
    std::size_t id(const std::string &name) const;
    const Material &material(std::size_t idx) const;
    const Material &material(const std::string &name) const;
    std::size_t size() const noexcept { return materials_.size(); }

private:
    std::vector<Material> materials_;
    std::unordered_map<std::string, std::size_t> lookup_;
};

} // namespace standalone
