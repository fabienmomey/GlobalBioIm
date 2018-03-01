% function m_opl_vmlmb_restore
%
%  ws = m_opl_vmlmb_restore(ws,x,f,grad);
%
%  Restore last line search starting point.  Calling this is only effective if                    
%    task is OPL_TASK_FG.
%
% See optimpacklegacy for explaination about VMLMB algorithm and its
% parameters
%
%	Definitions for optimization routines implemented in OptimPack
%	library.
%
%-----------------------------------------------------------------------------
%
%      Copyright (c) 2018, Ferreol SOULEZ.
%
%	This file is part of OptimPack.
%
%	OptimPack is  free software; you can redistribute  it and/or modify
%	it under the  terms of the GNU General  Public License as published
%	by the Free  Software Foundation; either version 2  of the License,
%	or (at your option) any later version.
%
%	OptimPack is  distributed in the hope  that it will  be useful, but
%	WITHOUT  ANY  WARRANTY;  without   even  the  implied  warranty  of
%	MERCHANTABILITY or  FITNESS FOR A PARTICULAR PURPOSE.   See the GNU
%	General Public License for more details.
%
%	You should have  received a copy of the  GNU General Public License
%	along with OptimPack (file  "LICENSE" in the top source directory);
%	if  not, write  to the  Free Software  Foundation, Inc.,  59 Temple
%	Place, Suite 330, Boston, MA 02111-1307 USA

