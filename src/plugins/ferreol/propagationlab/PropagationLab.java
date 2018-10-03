
/*  Copyright (C) 2017  Ferreol Soulez ferreol.soulez@univ-lyon1.fr
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *  */

package plugins.ferreol.propagationlab;

import icy.plugin.abstract_.Plugin;
import icy.plugin.interface_.PluginLibrary;

/**
 * Development PropagationLab
 *
 */
public class PropagationLab extends Plugin implements PluginLibrary {
    @SuppressWarnings("javadoc")
    public static void main(String [] args) {
        System.out.println("PropagationLab of collection of Fourier optic plugins\n"
                + "See https://github.com/FerreolS/PropagationLab for the  source   code\n"
                + "\n");
    }
}
