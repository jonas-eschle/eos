/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014 Stephan Jahn
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

#ifndef EOS_GUARD_SRC_UTILS_LIKELIHOOD_C_HH
#define EOS_GUARD_SRC_UTILS_LIKELIHOOD_C_HH 1

#include <eos/utils/log_likelihood.hh>

/// C like interface to handle a LogLikelihood object
extern "C" {
    using namespace eos;

    ///  constructor
    LogLikelihood * EOS_LogLikelihood_new();

    /// destructor
    void EOS_LogLikelihood_delete(LogLikelihood * ll);

    /*!
     * Add a constraint.
     *
     * @note The returned pointer must explicitly be deleted by using \c free.
     */
    char * EOS_LogLikelihood_add_constraint_by_name(LogLikelihood * ll,
                                                    const char * constraint_name,
                                                    const Options * options);
}

#endif
