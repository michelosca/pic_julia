# Copyright (C) 2021 Michel Osca Engelbrecht
#
# This file is part of PIC Julia.
#
# PIC Julia is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PIC Julia is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PIC Julia. If not, see <https://www.gnu.org/licenses/>.

module EvaluateExpressions

using Constants: e, amu, kb, epsilon_0, me
using SharedData: System, Waveform

function ReplaceExpressionValues(expr::Union{Int64, Float64, Expr},
    system::System; pos::Float64 = 0.0, waveform::Waveform=Waveform())

    copy_expr = copy(expr)
    if typeof(copy_expr)==Expr
        ReplaceSymbolValues!(copy_expr, system, pos, waveform)
    end

    value = eval(copy_expr)
    return value
end


function ReplaceSymbolValues!(expr::Expr, system::System, pos_x::Float64,
    waveform::Waveform)

    # system parameters
    ReplaceSystemValues!(expr, system)

    # Constant values
    ReplaceConstantValues!(expr)

    # Replace position value
    ReplaceSymbol!(expr, :x, pos_x)

    # Replace waveform value
    ReplaceWaveformValues!(expr, waveform)

end


function ReplaceSystemValues!(expr::Expr, system::System)
    # System parameters
    ReplaceSymbol!(expr, :Lx, system.Lx)
    ReplaceSymbol!(expr, :dx, system.dx)
    ReplaceSymbol!(expr, :x_min, system.x_min)
    ReplaceSymbol!(expr, :x_max, system.x_max)
    ReplaceSymbol!(expr, :time, system.time)
    ReplaceSymbol!(expr, :t_end, system.t_end)
    ReplaceSymbol!(expr, :dt, system.dt)
end


function ReplaceConstantValues!(expr::Expr)
    # Constant values
    ReplaceSymbol!(expr, :m_Ar, Expr(:call,:*,amu, 40 ))
    ReplaceSymbol!(expr, :amu, amu)
    ReplaceSymbol!(expr, :e, e)
    ReplaceSymbol!(expr, :kb, kb)
    ReplaceSymbol!(expr, :me, me)
    ReplaceSymbol!(expr, :eps0, epsilon_0)
end


function ReplaceSymbol!(expression::Expr, old, new)
    n = length(expression.args)
    for i in 1:n
        arg = expression.args[i]
        if (typeof(arg) == Expr)
            ReplaceSymbol!(arg, old, new)
        elseif (arg == old)
            expression.args[i] = new
        end
    end
end


function ReplaceWaveformValues!(expr::Expr, waveform::Waveform)
    # System parameters
    ReplaceSymbol!(expr, :f, waveform.freq)
end

end