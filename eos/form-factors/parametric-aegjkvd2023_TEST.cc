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

#include <test/test.hh>
#include <eos/form-factors/parametric-aegjkvd2023.hh>

#include <cmath>
#include <limits>
#include <vector>

using namespace test;
using namespace eos;

class ParametricAEGJKvD2023Test :
    public TestCase
{
    public:
        ParametricAEGJKvD2023Test() :
            TestCase("parametric_aegjkvd2023_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-7;

            Parameters p = Parameters::Defaults();
            p["mass::pi^+"]                    = 0.13957;
            p["pi->pi::t_0@AEGJKvD2023"]       = 0.0;
            p["pi->pi::a_(+,1)^0@AEGJKvD2023"] = 0.5;
            p["pi->pi::a_(+,1)^1@AEGJKvD2023"] = 0.25;
            p["pi->pi::a_(+,1)^2@AEGJKvD2023"] = 0.125;
            p["pi->pi::a_(+,1)^3@AEGJKvD2023"] = 0.0625;
            p["pi->pi::a_(+,1)^4@AEGJKvD2023"] = 0.03125;
            p["pi->pi::a_(+,1)^5@AEGJKvD2023"] = 0.015625;
            p["pi->pi::a_(+,1)^6@AEGJKvD2023"] = 7.8125e-3;
            p["pi->pi::a_(+,1)^7@AEGJKvD2023"] = 3.90625e-3;
            p["pi->pi::a_(+,1)^8@AEGJKvD2023"] = 1.953125e-3;
            p["pi->pi::a_(+,1)^9@AEGJKvD2023"] = 9.765625e-4;

            /* bounds */
            {
                AEGJKvD2023UnitarityBounds bounds(p, Options{});

                TEST_CHECK_NEARLY_EQUAL(bounds.bound_1m_I1(), 0.3333330, eps);
            }

            /* P->P factory */
            {
                std::shared_ptr<FormFactors<PToP>> ff = FormFactorFactory<PToP>::create("pi->pi::AEGJKvD2023", p, Options{ });

                TEST_CHECK(nullptr != ff);
            }

            /* f_+ and its auxilliary functions at spacelike and lightlike q2 <= 0.0 */
            {
                AEGJKvD2023FormFactors<PToP> ff(p, Options{ });

                const double chi = 3.52e-3;

                TEST_CHECK_NEARLY_EQUAL(ff.z(-1.0), 0.5762158, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.z( 0.0), 0.0,       eps);

                TEST_CHECK_NEARLY_EQUAL(ff.phi_p(ff.z(-1.0), chi), 1.30330542e-1, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.phi_p(ff.z( 0.0), chi), 2.96908774e-2, eps);

                TEST_CHECK_NEARLY_EQUAL(ff.f_p(-1.0),  2.47145192, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p( 0.0),  7.78799793, eps);
            }

            /* 0->PP factory */
            {
                std::shared_ptr<FormFactors<VacuumToPP>> ff = FormFactorFactory<VacuumToPP>::create("0->pipi::AEGJKvD2023", p, Options{ });

                TEST_CHECK(nullptr != ff);
            }

            /* f_+ at timlike q2 > 0.0 */
            {
                AEGJKvD2023FormFactors<VacuumToPP> ff(p, Options{ });

                const double chi = 3.52e-3;

                TEST_CHECK_NEARLY_EQUAL(real(ff.z( 0.0)),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(imag(ff.z( 0.0)),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(real(ff.z(+0.1)), -0.5583827, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(ff.z(+0.1)),  0.8295834,  eps);
                TEST_CHECK_NEARLY_EQUAL(real(ff.z(+0.5)),  0.6883234,  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(ff.z(+0.5)),  0.7254039,  eps);

                TEST_CHECK_NEARLY_EQUAL(real(ff.phi_p(ff.z( 0.0), chi)),  0.02969088, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(ff.phi_p(ff.z( 0.0), chi)),  0.0,        eps);
                TEST_CHECK_NEARLY_EQUAL(real(ff.phi_p(ff.z(+0.1), chi)),  0.00361219, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(ff.phi_p(ff.z(+0.1), chi)),  0.00827876, eps);
                TEST_CHECK_NEARLY_EQUAL(real(ff.phi_p(ff.z(+0.5), chi)), -0.0472899,  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(ff.phi_p(ff.z(+0.5), chi)),  0.0617105,  eps);

                TEST_CHECK_NEARLY_EQUAL(real(ff.f_p( 0.0)),   7.78799793, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p( 0.0)),   0.0,        eps);
                TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(+0.1)),  11.2844315,  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(+0.1)), -15.2218457,  eps);
                TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(+0.5)),  -0.4981299,  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(+0.5)),  -3.7979967,  eps);
            }
        }
} parametric_aegjkvd2023_test;
