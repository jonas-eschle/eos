/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2023 Danny van Dyk
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

#include <eos/b-decays/bq-to-dq-psd.hh>
#include <eos/form-factors/mesonic.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <string>

namespace eos
{
    using std::norm;

    /*
     * Decay: B_q -> D_q P, cf. [BBNS:2000A] (class I only, P = pi^- or K^-)
     */
    template <>
    struct Implementation<BqToDqPseudoscalar>
    {
        SpecifiedOption opt_model;

        std::shared_ptr<Model> model;

        UsedParameter hbar;

        UsedParameter g_fermi;

        QuarkFlavorOption opt_q;

        UsedParameter m_B;

        UsedParameter f_B;

        UsedParameter tau_B;

        UsedParameter m_D;

        UsedParameter f_D;

        UsedParameter m_P;

        UsedParameter f_P;

        SpecifiedOption opt_ff;

        std::shared_ptr<FormFactors<PToP>> ff;

        SpecifiedOption opt_cp_conjugate;

        bool cp_conjugate;

        UsedParameter mu;

        std::function<complex<double> ()> ckm_factor;
        std::function<WilsonCoefficients<bern::ClassIII> (bool)> wc;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            opt_model(o, options, "model"),
            model(Model::make(opt_model.value(), p, o)),
            opt_q(o, options, "q"),
            hbar(p["QM::hbar"], u),
            g_fermi(p["WET::G_Fermi"], u),
            m_B(p["mass::B_" + opt_q.str()], u),
            f_B(p["decay-constant::B_" + opt_q.str()], u),
            m_D(p["mass::D_" + opt_q.str()], u),
            f_D(p["decay-constant::D_" + opt_q.str()], u),
            m_P(p["mass::" + stringify(opt_q.value() == QuarkFlavor::down ? "K_u" : "pi^+")], u),
            f_P(p["decay-constant::" + stringify(opt_q.value() == QuarkFlavor::down ? "K_u" : "pi")], u),
            tau_B(p["life_time::B_" + opt_q.str()], u),
            opt_ff(o, options, "form-factors"),
            opt_cp_conjugate(o, options, "cp-conjugate"),
            cp_conjugate(destringify<bool>(opt_cp_conjugate.value())),
            mu(p[opt_q.str() + "bcu::mu"], u)
        {
            Context ctx("When constructing B_q->D_q P observable");

            switch (opt_q.value())
            {
                case QuarkFlavor::down:
                    ckm_factor = [this]() { return conj(model->ckm_ud()) * model->ckm_cb(); };
                    wc         = [this](const bool & cp_conjugate) { return model->wet_dbcu(cp_conjugate); };
                    ff         = FormFactorFactory<PToP>::create("B->D::" + opt_ff.value(), p, o);
                    break;
                case QuarkFlavor::strange:
                    ckm_factor = [this]() { return conj(model->ckm_us()) * model->ckm_cb(); };
                    wc         = [this](const bool & cp_conjugate) { return model->wet_sbcu(cp_conjugate); };
                    ff         = FormFactorFactory<PToP>::create("B_s->D_s::" + opt_ff.value(), p, o);
                    break;
                default:
                    throw InternalError("Invalid quark flavor: " + stringify(opt_q.value()));
            }
            u.uses(*model);
            u.uses(*ff);
        }

        complex<double> a1() const
        {
            const WilsonCoefficients<bern::ClassIII> wc = this->wc(cp_conjugate);

            return 1.0;
        }

        double decay_width() const
        {
            // cf. [BBNS:2000A], eq. (210), p. 80
            const complex<double> amplitude = g_fermi() / sqrt(2.0) * ckm_factor() * f_P() * ff->f_0(m_P * m_P)
                * (m_B * m_B - m_D * m_D) * this->a1();
            // cf. [BBNS:2000A], eq. (216), p. 80
            const double breakup_momentum = sqrt(lambda(m_B * m_B, m_D * m_D, m_P * m_P)) / (2.0 * m_B);

            // cf. [BBNS:2000A], eq. (221), p. 81
            return norm(amplitude) * breakup_momentum / (8.0 * M_PI);
        }

        double branching_ratio() const
        {
            return decay_width() * tau_B / hbar;
        }
    };

    const std::vector<OptionSpecification>
    Implementation<BqToDqPseudoscalar>::options
    {
        Model::option_specification(),
        FormFactorFactory<PToP>::option_specification(),
        { "cp-conjugate", { "true", "false" },  "false" },
        { "q",            { "s", "d" }                  },
        { "P",            { "K", "pi"}                  }
    };

    BqToDqPseudoscalar::BqToDqPseudoscalar(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BqToDqPseudoscalar>(new Implementation<BqToDqPseudoscalar>(parameters, options, *this))
    {
    }

    BqToDqPseudoscalar::~BqToDqPseudoscalar()
    {
    }

    double
    BqToDqPseudoscalar::branching_ratio() const
    {
        return _imp->branching_ratio();
    }

    double
    BqToDqPseudoscalar::decay_width() const
    {
        return _imp->decay_width();
    }

    const std::set<ReferenceName>
    BqToDqPseudoscalar::references
    {
        "BBNS:2000A"_rn
    };

    std::vector<OptionSpecification>::const_iterator
    BqToDqPseudoscalar::begin_options()
    {
        return Implementation<BqToDqPseudoscalar>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BqToDqPseudoscalar::end_options()
    {
        return Implementation<BqToDqPseudoscalar>::options.cend();
    }
}
