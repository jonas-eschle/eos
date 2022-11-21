/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2020-2022 Danny van Dyk
 * Copyright (c) 2020-2021 Eike S. Eberhard
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <eos/form-factors/parametric-aegjkvd2023.hh>
#include <eos/maths/szego-polynomial.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/options.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/qualified-name.hh>
#include <eos/utils/stringify.hh>

#include <numeric>

namespace eos
{
    const std::vector<OptionSpecification>
    AEGJKvD2023UnitarityBounds::_option_specifications;

    std::string
    AEGJKvD2023UnitarityBounds::_par_name(const std::string & ff, const std::string & isospin, const std::string & index) const
    {
        return "pi->pi::a_(" + ff + "," + isospin + ")^" + index + "@AEGJKvD2023";
    }

    AEGJKvD2023UnitarityBounds::AEGJKvD2023UnitarityBounds(const Parameters & p, const Options & /*o*/) :
        _a_fp_I1{{
            UsedParameter(p[_par_name("+", "1", "0")], *this),
            UsedParameter(p[_par_name("+", "1", "1")], *this),
            UsedParameter(p[_par_name("+", "1", "2")], *this),
            UsedParameter(p[_par_name("+", "1", "3")], *this),
            UsedParameter(p[_par_name("+", "1", "4")], *this),
            UsedParameter(p[_par_name("+", "1", "5")], *this),
            UsedParameter(p[_par_name("+", "1", "6")], *this),
            UsedParameter(p[_par_name("+", "1", "7")], *this),
            UsedParameter(p[_par_name("+", "1", "8")], *this),
            UsedParameter(p[_par_name("+", "1", "9")], *this)
        }}
    {
    }

    AEGJKvD2023UnitarityBounds::~AEGJKvD2023UnitarityBounds() = default;

    double
    AEGJKvD2023UnitarityBounds::bound_1m_I1() const
    {
        double sum = 0.0;
        for (auto & a : _a_fp_I1)
        {
            sum += a * a;
        }

        return sum;
    }

    const std::set<ReferenceName>
    AEGJKvD2023UnitarityBounds::references
    {
        "BL:1998A"_rn,
        "AEGJKvD:2023A"_rn
    };

    std::vector<OptionSpecification>::const_iterator
    AEGJKvD2023UnitarityBounds::begin_options()
    {
        return _option_specifications.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    AEGJKvD2023UnitarityBounds::end_options()
    {
        return _option_specifications.cend();
    }

    namespace aegjkvd2023
    {
        static const SzegoPolynomial<9u> szego_polynomial
        {
            64.0 / (5.0 * M_PI),
            {{
                +3.0 / 7.0,
                -5.0 / 9.0,
                +3.0 / 11.0,
                -5.0 / 13.0,
                +3.0 / 15.0,
                -5.0 / 17.0,
                +3.0 / 19.0,
                -5.0 / 21.0,
                +3.0 / 23.0
            }}
        };
    }

    AEGJKvD2023FormFactors<VacuumToPP>::AEGJKvD2023FormFactors(const Parameters & p, const Options & /*o*/) :
        _a_fp_I1{{
            UsedParameter(p[_par_name("+", "1", "0")], *this),
            UsedParameter(p[_par_name("+", "1", "1")], *this),
            UsedParameter(p[_par_name("+", "1", "2")], *this),
            UsedParameter(p[_par_name("+", "1", "3")], *this),
            UsedParameter(p[_par_name("+", "1", "4")], *this),
            UsedParameter(p[_par_name("+", "1", "5")], *this),
            UsedParameter(p[_par_name("+", "1", "6")], *this),
            UsedParameter(p[_par_name("+", "1", "7")], *this),
            UsedParameter(p[_par_name("+", "1", "8")], *this),
            UsedParameter(p[_par_name("+", "1", "9")], *this)
        }},
        _m_pi(p["mass::pi^+"], *this),
        _t_0(p["pi->pi::t_0@AEGJKvD2023"], *this)
    {
    }

    AEGJKvD2023FormFactors<VacuumToPP>::~AEGJKvD2023FormFactors() = default;

    FormFactors<VacuumToPP> *
    AEGJKvD2023FormFactors<VacuumToPP>::make(const Parameters & p, const Options & o)
    {
        return new AEGJKvD2023FormFactors<VacuumToPP>(p, o);
    }

    complex<double>
    AEGJKvD2023FormFactors<VacuumToPP>::z(const double & q2) const
    {
        return _z(q2, this->_t_0());
    }

    complex<double>
    AEGJKvD2023FormFactors<VacuumToPP>::phi_p(const complex<double> & z, const double & chi) const
    {
        const double t_p     = this->_t_p();
        const double t_0     = this->_t_0();
        const double tfactor = 1.0 - t_0 / t_p;
        const double Q2      = 2.0;

        // cf. [BL:1998A], eq. (4.6), p. 9.
        // note that the asymptotic factor ``(1.0 + z)^2 * (1.0 - z)^(-1/2)`` has been cancelled against the factor
        // in the redefined series expansion, i.e., ``(1.0 + z)^2 * (1.0 - z)^(-1/2) \sum_n b_n p_n(z) = \sum_n a_n z^n``.
        return 1.0 / sqrt(12.0 * M_PI * t_p * chi)
            * pow(tfactor, 5.0 / 4.0) * pow(sqrt(tfactor) * (1.0 + z) + (1.0 - z), -0.5)
            * pow(sqrt(1.0 + Q2 / t_p) * (1.0 - z) + sqrt(tfactor) * (1.0 + z), -3.0);
    }

    complex<double>
    AEGJKvD2023FormFactors<VacuumToPP>::series_p(const complex<double> & z, const std::array<double, 10u> & c) const
    {
        return std::inner_product(c.cbegin(), c.cend(), aegjkvd2023::szego_polynomial(z).cbegin(), complex<double>{ 0.0, 0.0 });
    }

    complex<double>
    AEGJKvD2023FormFactors<VacuumToPP>::f_p(const double & q2) const
    {
        // prepare expansion coefficients
        std::array<double, 10> a;
        std::copy(_a_fp_I1.cbegin(), _a_fp_I1.cend(), a.begin());

        const auto z           = this->z(q2);
        const auto chi         = 3.52e-3; // GeV^-2, cf. [BL:1998A], p. 13
        const auto phi         = this->phi_p(z, chi);
        const auto series      = this->series_p(z, a);

        return series / phi;
    }

    complex<double>
    AEGJKvD2023FormFactors<VacuumToPP>::f_t(const double & /*q2*/) const
    {
        throw InternalError("Not implemented!");
        return 0.0;
    }

    complex<double>
    AEGJKvD2023FormFactors<VacuumToPP>::f_0(const double & /*q2*/) const
    {
        return 0.0; // vanishes exactly
    }

    AEGJKvD2023FormFactors<PToP>::AEGJKvD2023FormFactors(const Parameters & p, const Options & /*o*/) :
        _a_fp_I1{{
            UsedParameter(p[_par_name("+", "1", "0")], *this),
            UsedParameter(p[_par_name("+", "1", "1")], *this),
            UsedParameter(p[_par_name("+", "1", "2")], *this),
            UsedParameter(p[_par_name("+", "1", "3")], *this),
            UsedParameter(p[_par_name("+", "1", "4")], *this),
            UsedParameter(p[_par_name("+", "1", "5")], *this),
            UsedParameter(p[_par_name("+", "1", "6")], *this),
            UsedParameter(p[_par_name("+", "1", "7")], *this),
            UsedParameter(p[_par_name("+", "1", "8")], *this),
            UsedParameter(p[_par_name("+", "1", "9")], *this)
        }},
        _m_pi(p["mass::pi^+"], *this),
        _t_0(p["pi->pi::t_0@AEGJKvD2023"], *this)
    {
    }

    AEGJKvD2023FormFactors<PToP>::~AEGJKvD2023FormFactors() = default;

    FormFactors<PToP> *
    AEGJKvD2023FormFactors<PToP>::make(const Parameters & p, const Options & o)
    {
        return new AEGJKvD2023FormFactors<PToP>(p, o);
    }

    double
    AEGJKvD2023FormFactors<PToP>::z(const double & q2) const
    {
        return _z(q2, this->_t_0());
    }

    double
    AEGJKvD2023FormFactors<PToP>::phi_p(const double & z, const double & chi) const
    {
        const double t_p     = this->_t_p();
        const double t_0     = this->_t_0();
        const double tfactor = 1.0 - t_0 / t_p;
        const double Q2      = 2.0;

        // cf. [BL:1998A], eq. (4.6), p. 9.
        // note that the asymptotic factor ``(1.0 + z)^2 * (1.0 - z)^(-1/2)`` has been cancelled against the factor
        // in the redefined series expansion, i.e., ``(1.0 + z)^2 * (1.0 - z)^(-1/2) \sum_n b_n p_n(z) = \sum_n a_n z^n``.
        return 1.0 / sqrt(12.0 * M_PI * t_p * chi)
            * pow(tfactor, 5.0 / 4.0) * pow(sqrt(tfactor) * (1.0 + z) + (1.0 - z), -0.5)
            * pow(sqrt(1.0 + Q2 / t_p) * (1.0 - z) + sqrt(tfactor) * (1.0 + z), -3.0);
    }

    double
    AEGJKvD2023FormFactors<PToP>::series_p(const double & z, const std::array<double, 10u> & c) const
    {
        return std::inner_product(c.cbegin(), c.cend(), aegjkvd2023::szego_polynomial(z).cbegin(), 0.0);
    }

    double
    AEGJKvD2023FormFactors<PToP>::f_p(const double & q2) const
    {
        // prepare expansion coefficients
        std::array<double, 10> a;
        std::copy(_a_fp_I1.cbegin(), _a_fp_I1.cend(), a.begin());

        const auto z           = this->z(q2);
        const auto chi         = 3.52e-3; // GeV^-2, cf. [BL:1998A], p. 13
        const auto phi         = this->phi_p(z, chi);
        // the asymptotics have been absorbed into the outer function to cancel a superficial divergence
        // as z -> +/- 1.0
        const auto series      = this->series_p(z, a);

        return series / phi;
    }

    double
    AEGJKvD2023FormFactors<PToP>::f_t(const double & /*q2*/) const
    {
        throw InternalError("Not implemented!");
        return 1.0;
    }

    double
    AEGJKvD2023FormFactors<PToP>::f_0(const double & /*q2*/) const
    {
        return 0.0; // vanishes exactly
    }

    template class AEGJKvD2023FormFactors<VacuumToPP>;
    template class AEGJKvD2023FormFactors<PToP>;
}
