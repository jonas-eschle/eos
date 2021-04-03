# vim: set sw=4 sts=4 et tw=120 :

# Copyright (c) 2018 Danny van Dyk
#
# This file is part of the EOS project. EOS is free software;
# you can redistribute it and/or modify it under the terms of the GNU General
# Public License version 2, as published by the Free Software Foundation.
#
# EOS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 59 Temple
# Place, Suite 330, Boston, MA  02111-1307  USA

import scipy

def __format_Parameter(p):
    name = ''
    latex = p.latex()
    if latex:
        name = r'$' + latex + r'$'
    else:
        name = p.name()

    return(
    """
    <table>
        <tr>
            <th>{name}</th>
            <td>(eos.Parameter)</td>
        </tr>
        <tr>
            <th>current value</th>
            <td><tt>{value}</tt></td>
        </tr>
        <tr>
            <th>default value</th>
            <td><tt>{central}</tt></td>
        </tr>
    </table>""".format(
        name=name,
        value=p.evaluate(),
        central=p.central()
    ))

def __format_KinematicVariable(kv):
    name = kv.name()

    return("""
        <table>
            <tr>
                <th>{name}</th>
                <td>(eos.KinematicVariable)</td>
            </tr>
            <tr>
                <th>current value</th>
                <td><tt>{value}</tt></td>
            </tr>
        </table>""".format(
            name=name,
            value=kv.evaluate()
    ))

def __format_Kinematics(k):
    result = '<table>\n'
    for kv in k:
        result += '<tr><th><tt style="color:grey">{name}</tt></th><td>{value}</td></tr>\n'.format(
            name=kv.name(),
            value=kv.evaluate()
        )
    result += '</table>\n'

    return result

def __format_Options(o):
    result = '<table>\n'
    for k, v in o:
        result += '<tr><th><tt style="color:grey">{key}</tt></th><td>{value}</td></tr>\n'.format(
            key=k,
            value=v
        )
    result += '</table>\n'

    return result

def __format_ObservableEntry(e):
    result = '<table>\n'
    result += '<tr><th>QualifedName</th><td><tt style="color:grey">{qn}</tt></td></tr>'.format(qn=e.name())
    result += '<tr><th>Description</th><td>$${latex}$$</td></tr>'.format(latex=e.latex())
    kvs = [kv for kv in e.kinematic_variables()]
    if len(kvs) > 0:
        result += '<tr><th rowspan={rows}>Kinematic Variables</th><td>{kv}</td></tr>'.format(rows=len(kvs),kv=kvs[0])
    for i in range(1, len(kvs)):
        result += '<tr><td>{kv}</td></tr>'.format(kv=kvs[i])
    result += '</table>'
    return(result)

def __format_Observable(obs):
    name = obs.name()
    kinematics = [(kv.name(), kv.evaluate()) for kv in obs.kinematics()]
    first_kinematics = "<td colspan=2>none</td>"
    further_kinematics = ""
    span_kinematics = 1
    if len(kinematics) > 0:
        first_kinematics = "<th>{kvn}</th><td>{kvv}</td>".format(kvn=kinematics[0][0], kvv=kinematics[0][1])
        further_kinematics = "\n".join([
            "<tr><th>{kvn}</th><td>{kvv}</td></tr>".format(kvn=kvn, kvv=kvv)
            for kvn, kvv in kinematics[1:]
        ])
        span_kinematics = len(kinematics)
    options = [(ok, ov) for ok, ov in obs.options()]
    first_options = "<td colspan=2>none</td>"
    further_options = ""
    span_options = 1
    if len(options) > 0:
        first_options = "<th>{ok}</th><td>{ov}</td>".format(ok=options[0][0], ov=options[0][1])
        further_options = "\n".join([
            "<tr><th>{ok}</th><td>{ov}</td></tr>".format(ok=ok, ov=ov)
            for ok, ov in options[1:]
        ])
        span_options = len(options)

    return("""
        <table>
            <tr>
                <th>{name}</th>
                <td colspan="2">(eos.Observable)</td>
            </tr>
            <tr>
                <th rowspan="{span_kinematics}">kinematics</th>
                {first_kinematics}
            </tr>
            {further_kinematics}
            <tr>
                <th rowspan="{span_options}">options</th>
                {first_options}
            </tr>
            {further_options}
            <tr>
                <th>current value</th>
                <td colspan="2"><tt>{value:1.4g}</tt></td>
            </tr>
        </table>""".format(
            name=name,
            value=obs.evaluate(),
            first_kinematics=first_kinematics,
            span_kinematics=span_kinematics,
            further_kinematics=further_kinematics,
            first_options=first_options,
            span_options=span_options,
            further_options=further_options
    ))

def __format_GoodnessOfFit(gof):
    result = '<table>\n'
    result += '<tr><th>constraint</th><th>&chi;<sup>2</sup></th><th>d.o.f.</th></tr>\n'
    for entry in gof:
        result += '<tr><td><tt>{name}</tt></td><td>{chi2:6.4f}</td><td>{dof}</td></tr>\n'.format(
            name=entry[0], chi2=entry[1].chi2, dof=entry[1].dof)
    result += '</table><br/>\n'
    chi2 = gof.total_chi_square()
    dof  = gof.total_degrees_of_freedom()
    pvalue = 1.0 - scipy.stats.chi2(dof).cdf(chi2)
    result += '<table>\n'
    result += '<tr><th>total &chi;<sup>2</sup></th><td>{chi2:6.4f}</td></tr>\n'.format(chi2=chi2)
    result += '<tr><th>total degrees of freedom</th><td>{dof}</td></tr>\n'.format(dof=dof)
    result += '<tr><th>p-value</th><td>{p:6.4f}%</td></tr>\n'.format(p=pvalue * 100)
    result += '</table>\n'
    return(result)
