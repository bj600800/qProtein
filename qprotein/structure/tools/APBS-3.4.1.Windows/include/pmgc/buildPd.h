/**
 *  @ingroup PMGC
 *  @author  Tucker Beck [fortran ->c translation], Michael Holst [original]
 *  @brief Builds prolongation matrix
 *  @version $Id:
 *
 *  @attention
 *  @verbatim
 *
 * APBS -- Adaptive Poisson-Boltzmann Solver
 *
 * Nathan A. Baker (nathan.baker@pnl.gov)
 * Pacific Northwest National Laboratory
 *
 * Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2010-2020 Battelle Memorial Institute. Developed at the Pacific Northwest National Laboratory, operated by Battelle Memorial Institute, Pacific Northwest Division for the U.S. Department Energy.  Portions Copyright (c) 2002-2010, Washington University in St. Louis.  Portions Copyright (c) 2002-2010, Nathan A. Baker.  Portions Copyright (c) 1999-2002, The Regents of the University of California. Portions Copyright (c) 1995, Michael Holst.
 * All rights reserved.
 *
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * -  Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * - Neither the name of Washington University in St. Louis nor the names of its
 * contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @endverbatim
 */

#ifndef _BUILDPD_H_
#define _BUILDPD_H_

#include "apbscfg.h"

#include "maloc/maloc.h"

#include "generic/vhal.h"
#include "generic/vmatrix.h"

VEXTERNC void VbuildP(
        int    *nxf,    ///< @todo: doc
        int    *nyf,    ///< @todo: doc
        int    *nzf,    ///< @todo: doc
        int    *nxc,    ///< @todo: doc
        int    *nyc,    ///< @todo: doc
        int    *nzc,    ///< @todo: doc
        int    *mgprol, ///< @todo: doc
        int    *ipc,    ///< @todo: doc
        double *rpc,    ///< @todo: doc
        double *pc,     ///< @todo: doc
        double *ac,     ///< @todo: doc
        double *xf,     ///< @todo: doc
        double *yf,     ///< @todo: doc
        double *zf      ///< @todo: doc
);

VEXTERNC void VbuildP_trilin(
        int *nxf,   ///< @todo: doc
        int *nyf,   ///< @todo: doc
        int *nzf,   ///< @todo: doc
        int *nxc,   ///< @todo: doc
        int *nyc,   ///< @todo: doc
        int *nzc,   ///< @todo: doc
        double *pc, ///< @todo: doc
        double *xf, ///< @todo: doc
        double *yf, ///< @todo: doc
        double *zf  ///< @todo: doc
);

VEXTERNC void VbuildPb_trilin(
        int    *nxf,  ///< @todo: doc
        int    *nyf,  ///< @todo: doc
        int    *nzf,  ///< @todo: doc
        int    *nxc,  ///< @todo: doc
        int    *nyc,  ///< @todo: doc
        int    *nzc,  ///< @todo: doc
        double *oPC,  ///< @todo: doc
        double *oPN,  ///< @todo: doc
        double *oPS,  ///< @todo: doc
        double *oPE,  ///< @todo: doc
        double *oPW,  ///< @todo: doc
        double *oPNE, ///< @todo: doc
        double *oPNW, ///< @todo: doc
        double *oPSE, ///< @todo: doc
        double *oPSW, ///< @todo: doc
        double *uPC,  ///< @todo: doc
        double *uPN,  ///< @todo: doc
        double *uPS,  ///< @todo: doc
        double *uPE,  ///< @todo: doc
        double *uPW,  ///< @todo: doc
        double *uPNE, ///< @todo: doc
        double *uPNW, ///< @todo: doc
        double *uPSE, ///< @todo: doc
        double *uPSW, ///< @todo: doc
        double *dPC,  ///< @todo: doc
        double *dPN,  ///< @todo: doc
        double *dPS,  ///< @todo: doc
        double *dPE,  ///< @todo: doc
        double *dPW,  ///< @todo: doc
        double *dPNE, ///< @todo: doc
        double *dPNW, ///< @todo: doc
        double *dPSE, ///< @todo: doc
        double *dPSW, ///< @todo: doc
        double *xf,   ///< @todo: doc
        double *yf,   ///< @todo: doc
        double *zf    ///< @todo: doc
);

VEXTERNC void VbuildP_op7(
        int    *nxf, ///< @todo: doc
        int    *nyf, ///< @todo: doc
        int    *nzf, ///< @todo: doc
        int    *nxc, ///< @todo: doc
        int    *nyc, ///< @todo: doc
        int    *nzc, ///< @todo: doc
        int    *ipc, ///< @todo: doc
        double *rpc, ///< @todo: doc
        double  *ac, ///< @todo: doc
        double  *pc  ///< @todo: doc
);

VEXTERNC void VbuildPb_op7(
                int     *nxf, ///< @todo: doc
                int     *nyf, ///< @todo: doc
                int     *nzf, ///< @todo: doc
                int     *nxc, ///< @todo: doc
                int     *nyc, ///< @todo: doc
                int     *nzc, ///< @todo: doc
                int     *ipc, ///< @todo: doc
                double  *rpc, ///< @todo: doc
                double   *oC, ///< @todo: doc
                double   *oE, ///< @todo: doc
                double   *oN, ///< @todo: doc
                double   *uC, ///< @todo: doc
                double  *oPC, ///< @todo: doc
                double  *oPN, ///< @todo: doc
                double  *oPS, ///< @todo: doc
                double  *oPE, ///< @todo: doc
                double  *oPW, ///< @todo: doc
                double *oPNE, ///< @todo: doc
                double *oPNW, ///< @todo: doc
                double *oPSE, ///< @todo: doc
                double *oPSW, ///< @todo: doc
                double  *uPC, ///< @todo: doc
                double  *uPN, ///< @todo: doc
                double  *uPS, ///< @todo: doc
                double  *uPE, ///< @todo: doc
                double  *uPW, ///< @todo: doc
                double *uPNE, ///< @todo: doc
                double *uPNW, ///< @todo: doc
                double *uPSE, ///< @todo: doc
                double *uPSW, ///< @todo: doc
                double  *dPC, ///< @todo: doc
                double  *dPN, ///< @todo: doc
                double  *dPS, ///< @todo: doc
                double  *dPE, ///< @todo: doc
                double  *dPW, ///< @todo: doc
                double *dPNE, ///< @todo: doc
                double *dPNW, ///< @todo: doc
                double *dPSE, ///< @todo: doc
                double *dPSW  ///< @todo: doc
);

VEXTERNC void VbuildP_op27(
                int      *nxf, ///< @todo: doc
                int      *nyf, ///< @todo: doc
                int      *nzf, ///< @todo: doc
                int      *nxc, ///< @todo: doc
                int      *nyc, ///< @todo: doc
                int      *nzc, ///< @todo: doc
                int      *ipc, ///< @todo: doc
                double   *rpc, ///< @todo: doc
                double    *ac, ///< @todo: doc
                double    *pc  ///< @todo: doc
);

VEXTERNC void VbuildPb_op27(
                int      *nxf, ///< @todo: doc
                int      *nyf, ///< @todo: doc
                int      *nzf, ///< @todo: doc
                int      *nxc, ///< @todo: doc
                int      *nyc, ///< @todo: doc
                int      *nzc, ///< @todo: doc
                int      *ipc, ///< @todo: doc
                double   *rpc, ///< @todo: doc
                double    *oC, ///< @todo: doc
                double    *oE, ///< @todo: doc
                double    *oN, ///< @todo: doc
                double    *uC, ///< @todo: doc
                double   *oNE, ///< @todo: doc
                double   *oNW, ///< @todo: doc
                double    *uE, ///< @todo: doc
                double    *uW, ///< @todo: doc
                double    *uN, ///< @todo: doc
                double    *uS, ///< @todo: doc
                double   *uNE, ///< @todo: doc
                double   *uNW, ///< @todo: doc
                double   *uSE, ///< @todo: doc
                double   *uSW, ///< @todo: doc
                double   *oPC, ///< @todo: doc
                double   *oPN, ///< @todo: doc
                double   *oPS, ///< @todo: doc
                double   *oPE, ///< @todo: doc
                double   *oPW, ///< @todo: doc
                double  *oPNE, ///< @todo: doc
                double  *oPNW, ///< @todo: doc
                double  *oPSE, ///< @todo: doc
                double  *oPSW, ///< @todo: doc
                double   *uPC, ///< @todo: doc
                double   *uPN, ///< @todo: doc
                double   *uPS, ///< @todo: doc
                double   *uPE, ///< @todo: doc
                double   *uPW, ///< @todo: doc
                double  *uPNE, ///< @todo: doc
                double  *uPNW, ///< @todo: doc
                double  *uPSE, ///< @todo: doc
                double  *uPSW, ///< @todo: doc
                double   *dPC, ///< @todo: doc
                double   *dPN, ///< @todo: doc
                double   *dPS, ///< @todo: doc
                double   *dPE, ///< @todo: doc
                double   *dPW, ///< @todo: doc
                double  *dPNE, ///< @todo: doc
                double  *dPNW, ///< @todo: doc
                double  *dPSE, ///< @todo: doc
                double  *dPSW  ///< @todo: doc
);

#endif /* _BUILDPD_H_ */
