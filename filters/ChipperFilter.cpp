/******************************************************************************
 * Copyright (c) 2010, Andrew Bell
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following
 * conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in
 *       the documentation and/or other materials provided
 *       with the distribution.
 *     * Neither the name of the Andrew Bell or libLAS nor the names of
 *       its contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGE.
 ****************************************************************************/

#include "ChipperFilter.hpp"

#include <iostream>
#include <limits>

/**
The objective is to split the region into non-overlapping blocks, each
containing approximately the same number of points, as specified by the
user.  We'd also like the blocks closer to squares (or cubes) than not.

First, the points are read into arrays - one for the x direction, one for
the y direction, and one for the z direction when use_z is set.  The arrays
are sorted and are initialized with indices into the other arrays of the
location of the other coordinate of the same point.

Partitions are created that place the maximum number of points in a
block, subject to the user-defined threshold, using a cumulate and round
procedure.

The distance of the point-space is checked in each direction and the
wider dimension is chosen for splitting at an appropriate partition point.
The points in the narrower directions are copied to locations in the spare
arrays at one side or the other of the chosen partition, and that portion
of the spare arrays then become the active arrays for the narrow directions.
This avoids resorting of the arrays, which are already sorted.

This procedure is then recursively applied to the created blocks until
they contain only one or two partitions.  In the case of one or two
partitions we are done, and we simply store away the contents of the
blocks.
**/

#include <pdal/pdal_macros.hpp>
#include <pdal/util/ProgramArgs.hpp>

namespace pdal
{

static PluginInfo const s_info = PluginInfo(
    "filters.chipper",
    "Organize points into spatially contiguous, squarish, and non-overlapping chips.",
    "http://pdal.io/stages/filters.chipper.html" );

CREATE_STATIC_PLUGIN(1, 0, ChipperFilter, Filter, s_info)

std::string ChipperFilter::getName() const { return s_info.name; }

void ChipperFilter::addArgs(ProgramArgs& args)
{
    args.add("capacity", "Maximum number of points per cell", m_threshold,
        (PointId) 5000u);
    args.add("use_z", "Whether Z values are considered", m_useZ, false);
}


PointViewSet ChipperFilter::run(PointViewPtr view)
{
    point_count_t point_count = view->size();
    if (point_count == 0)
        return m_outViews;

    m_inView = view;
    load(*view.get(), m_xvec, m_yvec, m_zvec, m_spare1, m_spare2);
    partition(point_count);
    if (m_zvec.size() == 0)
        decideSplit2D(m_xvec, m_yvec, m_spare1, 0, m_partitions.size() - 1);
    else
        decideSplit3D(m_xvec, m_yvec, m_zvec, m_spare1, m_spare2, 0, m_partitions.size() - 1);
    return m_outViews;
}


void ChipperFilter::load(PointView& view, ChipRefList& xvec, ChipRefList& yvec,
    ChipRefList& zvec, ChipRefList& spare1, ChipRefList& spare2)
{
    point_count_t idx;
    std::vector<ChipPtRef>::iterator it;

    xvec.reserve(view.size());
    yvec.reserve(view.size());
    spare1.resize(view.size());
    if (m_useZ) {
        zvec.reserve(view.size());
        spare2.resize(view.size());
    }

    for (PointId i = 0; i < view.size(); ++i)
    {
        ChipPtRef xref;

        xref.m_pos = view.getFieldAs<double>(Dimension::Id::X, i);
        xref.m_ptindex = i;
        xvec.push_back(xref);

        ChipPtRef yref;

        yref.m_pos = view.getFieldAs<double>(Dimension::Id::Y, i);
        yref.m_ptindex = i;
        yvec.push_back(yref);

        ChipPtRef zref;

        if (m_useZ) {
            zref.m_pos = view.getFieldAs<double>(Dimension::Id::Z, i);
            zref.m_ptindex = i;
            zvec.push_back(zref);
        }
    }

    // Sort xvec
    std::stable_sort(xvec.begin(), xvec.end());

    // Set the other index in yvec and zvec
    for (size_t i = 0; i < xvec.size(); ++i)
    {
        idx = xvec[i].m_ptindex;
        if (m_useZ) {
            yvec[idx].m_o2index = i;
            zvec[idx].m_o1index = i;
        } else {
            yvec[idx].m_o1index = i;
        }
    }

    // Sort yvec
    std::stable_sort(yvec.begin(), yvec.end());

    // Set the other index in xvec and zvec
    for (size_t i = 0; i < yvec.size(); ++i)
    {
        if (m_useZ) {
            xvec[yvec[i].m_o2index].m_o1index = i;
            idx = yvec[i].m_ptindex;
            zvec[idx].m_o2index = i;
        } else {
            xvec[yvec[i].m_o1index].m_o1index = i;
        }
    }

    // Sort zvec
    std::stable_sort(zvec.begin(), zvec.end());

    // Set the other index in xvec and yvec
    for (size_t i = 0; i < zvec.size(); ++i)
    {
        xvec[zvec[i].m_o1index].m_o2index = i;
        yvec[zvec[i].m_o2index].m_o1index = i;
    }
}




#ifdef _WIN32
inline long lround(double d)
{
	long l;

	if (d < 0)
		l = (long)ceil(d - .5);
	else
		l = (long)floor(d + .5);
	return l;
}
#endif


// Build a list of partitions.  The partition is the size of each block in
// the x, y (and z) directions in number of points.
void ChipperFilter::partition(point_count_t size)
{
    size_t num_partitions;

    num_partitions = size / m_threshold;
    if (size % m_threshold)
        num_partitions++;

    // This is a standard statistics cumulate and round.  It distributes
    // the points into partitions such the "extra" points are reasonably
    // distributed among the partitions.
    double total(0.0);
    double partition_size = static_cast<double>(size) / num_partitions;
    m_partitions.push_back(0);
    for (size_t i = 0; i < num_partitions; ++i)
    {
        total += partition_size;
        size_t itotal = lround(total);
        m_partitions.push_back(itotal);
    }
}

void ChipperFilter::decideSplit2D(ChipRefList& v1, ChipRefList& v2,
        ChipRefList& spare, PointId pleft, PointId pright)
{
    double v1range;
    double v2range;
    uint32_t left = m_partitions[pleft];
    uint32_t right = m_partitions[pright] - 1;

    // Decide the wider direction of the block, and split in that direction
    // to maintain squareness.
    v1range = v1[right].m_pos - v1[left].m_pos;
    v2range = v2[right].m_pos - v2[left].m_pos;

    if (v1range > v2range)
        split2D(v1, v2, spare, pleft, pright);
    else
        split2D(v2, v1, spare, pleft, pright);
}

void ChipperFilter::split2D(ChipRefList& wide, ChipRefList& narrow,
        ChipRefList& spare, PointId pleft, PointId pright)
{
    PointId lstart;
    PointId rstart;
    PointId pcenter;
    PointId left;
    PointId right;
    PointId center;

    left = m_partitions[pleft];
    right = m_partitions[pright] - 1;

    // There are two cases in which we are done.
    // 1) We have a distance of one between left and right.
    // 2) We have a distance of two between left and right.

    if (pright - pleft == 1)
        emit(wide, left, right);
    else if (pright - pleft == 2) {
        center = m_partitions[pright - 1];
        emit(wide,
             left,
             center - 1);
        emit(wide,
             center,
             right);
    } else {
        pcenter = (pleft + pright) / 2;
        center = m_partitions[pcenter];

        // We are splitting in the wide direction - split elements in the
        // narrow array by copying them to the spare array in the correct
        // partition.  The spare array then becomes the active narrow array
        // for the [left,right] partition.
        lstart = left;
        rstart = center;
        for (PointId i = left; i <= right; ++i)
        {
            if (narrow[i].m_o1index < center)
            {
                spare[lstart] = narrow[i];
                wide[narrow[i].m_o1index].m_o1index = lstart;
                lstart++;
            }
            else
            {
                spare[rstart] = narrow[i];
                wide[narrow[i].m_o1index].m_o1index = rstart;
                rstart++;
            }
        }

        decideSplit2D(wide, spare, narrow, pleft, pcenter);
        decideSplit2D(wide, spare, narrow, pcenter, pright);
    }
}

void ChipperFilter::decideSplit3D(ChipRefList& v1, ChipRefList& v2, ChipRefList& v3,
    ChipRefList& spare1, ChipRefList& spare2, PointId pleft, PointId pright)
{
    double v1range;
    double v2range;
    double v3range;
    uint32_t left = m_partitions[pleft];
    uint32_t right = m_partitions[pright] - 1;

    // Decide the wider direction of the block, and split in that direction
    // to maintain cubeness.
    v1range = v1[right].m_pos - v1[left].m_pos;
    v2range = v2[right].m_pos - v2[left].m_pos;
    v3range = v3[right].m_pos - v3[left].m_pos;

    if (v1range > v2range) {
        if (v1range > v3range) {
            split3D(v1, v2, v3, spare1, spare2, pleft, pright);
        } else {
            split3D(v3, v1, v2, spare1, spare2, pleft, pright);
        }
    } else {
        if (v2range > v3range) {
            split3D(v2, v3, v1, spare1, spare2, pleft, pright);
        } else {
            split3D(v3, v1, v2, spare1, spare2, pleft, pright);
        }
    }
}

void ChipperFilter::split3D(ChipRefList& wide, ChipRefList& narrow1, ChipRefList& narrow2,
    ChipRefList& spare1, ChipRefList& spare2, PointId pleft, PointId pright)
{
    PointId lstart;
    PointId rstart;
    PointId pcenter;
    PointId left;
    PointId right;
    PointId center;

    left = m_partitions[pleft];
    right = m_partitions[pright] - 1;

    // There are two cases in which we are done.
    // 1) We have a distance of one between left and right.
    // 2) We have a distance of two between left and right.

    if (pright - pleft == 1)
        emit(wide, left, right);
    else if (pright - pleft == 2) {
        center = m_partitions[pright - 1];
        emit(wide,
             left,
             center - 1);
        emit(wide,
             center,
             right);
    } else {
        pcenter = (pleft + pright) / 2;
        center = m_partitions[pcenter];


        // We are splitting in the wide direction - split elements in the
        // narrow arrays by copying them to the spare arrays in the correct
        // partition.  The spare arrays then become the active narrow arrays
        // for the [left,right] partition.
        lstart = left;
        rstart = center;
        for (PointId i = left; i <= right; ++i)
        {
            if (narrow1[i].m_o2index < center)
            {
                spare1[lstart] = narrow1[i];
                wide[narrow1[i].m_o2index].m_o1index = lstart;
                narrow2[narrow1[i].m_o1index].m_o2index = lstart;
                lstart++;
            }
            else
            {
                spare1[rstart] = narrow1[i];
                wide[narrow1[i].m_o2index].m_o1index = rstart;
                narrow2[narrow1[i].m_o1index].m_o2index = rstart;
                rstart++;
            }
        }

        lstart = left;
        rstart = center;
        for (PointId i = left; i <= right; ++i) {
            if (narrow2[i].m_o1index < center)
            {
                spare2[lstart] = narrow2[i];
                wide[narrow2[i].m_o1index].m_o2index = lstart;
                spare1[narrow2[i].m_o2index].m_o1index = lstart;
                lstart++;
            }
            else
            {
                spare2[rstart] = narrow2[i];
                wide[narrow2[i].m_o1index].m_o2index = rstart;
                spare1[narrow2[i].m_o2index].m_o1index = rstart;
                rstart++;
            }
        }

        decideSplit3D(wide, spare1, spare2, narrow1, narrow2, pleft, pcenter);
        decideSplit3D(wide, spare1, spare2, narrow1, narrow2, pcenter, pright);
    }
}

void ChipperFilter::emit(ChipRefList& wide, PointId widemin, PointId widemax)
{
    PointViewPtr view = m_inView->makeNew();
    for (size_t idx = widemin; idx <= widemax; ++idx)
        view->appendPoint(*m_inView.get(), wide[idx].m_ptindex);

    m_outViews.insert(view);
}

} // namespace pdal

