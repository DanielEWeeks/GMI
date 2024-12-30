#   GMI: Genetic Map Interpolator
#   Copyright (C) 2010 Nandita Mukhopadhyay, Xinyu Tang, Daniel E. Weeks
# 
#   This file is part of the GMI program, which is free software; you
#   can redistribute it and/or modify it under the terms of the GNU
#   General Public License as published by the Free Software Foundation;
#   either version 2 of the License, or (at your option) any later
#   version.
# 
#   GMI is distributed in the hope that it will be useful, but WITHOUT
#   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#   for more details.
# 
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# 
#   For further information contact:
#       Daniel E. Weeks
#       e-mail: weeks@pitt.edu
#
# ===========================================================================

use IO::File;
use File::Copy;

# Declare the subroutines

# Perl trim function to remove whitespace from the start and end of the string

sub trim($)
{
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

sub check_installed_utils {
    
    my $prog = shift;
    my $result = `which $prog`;

    if ($result) {
	return 1;
    }
    else {
	return 0;
    }
}

1;
