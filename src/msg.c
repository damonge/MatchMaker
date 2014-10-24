///////////////////////////////////////////////////////////////////////
//                                                                   //
//   Copyright 2012 David Alonso                                     //
//                                                                   //
//                                                                   //
// This file is part of MatchMaker.                                  //
//                                                                   //
// MatchMaker is free software: you can redistribute it and/or       //
// modify it under the terms of the GNU General Public License as    //
// published by the Free Software Foundation, either version 3 of    //
// the License, or (at your option) any later version.               //
//                                                                   //
// MatchMaker is distributed in the hope that it will be useful, but //
// WITHOUT ANY WARRANTY; without even the implied warranty of        //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU //
// General Public License for more details.                          //
//                                                                   //
// You should have received a copy of the GNU General Public License //
// along with MatchMaker.  If not, see                               //
//        <http://www.gnu.org/licenses/>.                            //
//                                                                   //
///////////////////////////////////////////////////////////////////////
//
// Utility functions to write message, record computation time ...
//

#include <stdio.h>
#include <stdarg.h>
#include <mpi.h>

#include "common.h"

static int myrank=-1;

// Initialize using msg_ functions.
void msg_init()
{
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
}

void msg_printf(const char *fmt, ...)
{
  if(myrank==0) {
    va_list argp;

    va_start(argp,fmt);
    vfprintf(stdout,fmt,argp);
    fflush(stdout);
    va_end(argp);
  }
}

void msg_abort(const int errret,const char *fmt, ...)
{
  va_list argp;

  va_start(argp,fmt);
  vfprintf(stderr,fmt,argp);
  va_end(argp);

  MPI_Abort(MPI_COMM_WORLD,errret);
}  
