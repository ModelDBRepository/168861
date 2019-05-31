#!/bin/bash
#
# File author(s):  <Julian Andres Garcia Grajales>, (C) 2014
#
# Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
#
#

mencoder mf://GUI_captures/*.bmp -mf fps=10 -o movie.avi -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=800
