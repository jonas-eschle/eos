/* vim: set sw=4 sts=4 et tw=120 foldmethod=syntax : */

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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_AEGJKVD2023_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_AEGJKVD2023_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/maths/complex.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/reference-name.hh>

#include <array>

namespace eos
{
    template <typename Transition_> class AEGJKvD2023FormFactors;

    template <> class AEGJKvD2023FormFactors<VacuumToPP> :
        public FormFactors<VacuumToPP>
    {
        private:
            // parameters for form factors f_+
            std::array<UsedParameter, 10u> _a_fp_I1;

            // hadron masses
            UsedParameter _m_pi;

            // parameter for zero point of z
            UsedParameter _t_0;

            inline std::string _par_name(const std::string & ff, const std::string & isospin, const std::string & index) const
            {
                return "pi->pi::a_(" + ff + "," + isospin + ")^" + index + "@AEGJKvD2023";
            }

            inline double _t_p() const
            {
                return 4.0 * _m_pi() * _m_pi();
            }

            inline complex<double> _z(const double & q2, const complex<double> & t_0) const
            {
                const complex<double> t_p = complex<double>{ _t_p(), 0.0};

                return (sqrt(t_p - q2) - sqrt(t_p - t_0)) / (sqrt(t_p - q2) + sqrt(t_p - t_0));
            }

        public:
            AEGJKvD2023FormFactors(const Parameters & p, const Options & o);
            ~AEGJKvD2023FormFactors();

            static FormFactors<VacuumToPP> * make(const Parameters & p, const Options & o);

            /* f_+ auxilliary functions */
            complex<double> z(const double & q2) const;
            complex<double> phi_p(const complex<double> & z, const double & chi) const;
            complex<double> series_p(const complex<double> & z, const std::array<double, 10u> & c) const;

            virtual complex<double> f_p(const double & q2) const;
            virtual complex<double> f_t(const double & q2) const;
            virtual complex<double> f_0(const double & q2) const;
    };

    extern template class AEGJKvD2023FormFactors<VacuumToPP>;

    template <> class AEGJKvD2023FormFactors<PToP> :
        public FormFactors<PToP>
    {
        private:
            // parameters for form factors f_+
            std::array<UsedParameter, 10u> _a_fp_I1;

            // hadron masses
            UsedParameter _m_pi;

            // parameter for zero point of z
            UsedParameter _t_0;

            inline std::string _par_name(const std::string & ff, const std::string & isospin, const std::string & index) const
            {
                return "pi->pi::a_(" + ff + "," + isospin + ")^" + index + "@AEGJKvD2023";
            }

            inline double _t_p() const
            {
                return 4.0 * _m_pi() * _m_pi();
            }

            inline double _z(const double & q2, const double & t_0) const
            {
                const auto t_p = _t_p();
                return (sqrt(t_p - q2) - sqrt(t_p - t_0)) / (sqrt(t_p - q2) + sqrt(t_p - t_0));
            }

        public:
            AEGJKvD2023FormFactors(const Parameters & p, const Options & o);
            ~AEGJKvD2023FormFactors();

            static FormFactors<PToP> * make(const Parameters & p, const Options & o);

            /* f_+ auxilliary functions */
            double z(const double & q2) const;
            double phi_p(const double & z, const double & chi) const;
            double series_p(const double & z, const std::array<double, 10u> & c) const;

            virtual double f_p(const double & q2) const override;
            virtual double f_t(const double & q2) const override;
            virtual double f_0(const double & q2) const override;
    };

    extern template class AEGJKvD2023FormFactors<PToP>;

    class AEGJKvD2023UnitarityBounds :
        public ParameterUser
    {
        private:
            static const std::vector<OptionSpecification> _option_specifications;

            // parameters for form factor f_+ with I = 1
            std::array<UsedParameter, 10u> _a_fp_I1;

            std::string _par_name(const std::string & ff, const std::string & isospin, const std::string & index) const;

        public:
            AEGJKvD2023UnitarityBounds(const Parameters & p, const Options & o);

            ~AEGJKvD2023UnitarityBounds();

            double bound_1m_I1() const;

            /*!
             * References used in the computation of our observables.
             */
            static const std::set<ReferenceName> references;

            /*!
             * Options used in the computation of our observables.
             */
            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();
    };
}

#endif
