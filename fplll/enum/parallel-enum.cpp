/* Copyright (C) 2019 Marc Stevens

   This file is part of fplll. fplll is free software: you
   can redistribute it and/or modify it under the terms of the GNU Lesser
   General Public License as published by the Free Software Foundation,
   either version 2.1 of the License, or (at your option) any later version.

   fplll is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with fplll. If not, see <http://www.gnu.org/licenses/>. */

#include "parallel-enum.h"

FPLLL_BEGIN_NAMESPACE

template <typename ZT, typename FT>
void ParallelEnumerationDyn<ZT, FT>::enumerate(int first, int last, FT &fmaxdist, long fmaxdistexpo,
                                               int split, const vector<FT> &target_coord,
                                               const vector<enumxt> &subtree,
                                               const vector<enumf> &pruning)
{
  if (last == -1)
    last = _gso.d;
  int d = last - first;
  if (split == -1)
    split = 1 + (d - subtree.size()) / 20;
  if (split < subtree.size() + 1 || split >= d - 1)
    throw std::runtime_error("ParallelEnumerationDyn::enumerate(): split out of bounds");

  /* prepare bottom enumeration jobs */
  _maxdist = -1;

  _first        = first;
  _last         = last;
  _split        = split;
  _d            = d;
  _fmaxdist     = fmaxdist;
  _fmaxdistexpo = fmaxdistexpo;
  _target_coord = target_coord;
  _subtree      = subtree;
  _pruning      = pruning;
  _dual         = dual;

  _toptrees.clear();
  _toptrees.reserve(1 << 16);

  _bottom_eval.clear();
  for (unsigned i = 0; i < get_threads(); ++i)
    _bottom_eval.emplace_back(*this);

  _bottom_enums.clear();
  for (unsigned i = 0; i < get_threads(); ++i)
    _bottom_enums.emplace_back(_gso, _bottom_eval[i]);

  _bottom_fmaxdist.clear();
  for (unsigned i = 0; i < get_threads(); ++i)
    _bottom_fmaxdist.emplace_back(fmaxdist);

  for (unsigned i = 0; i < get_threads() - 1; ++i)
    threadpool.push([this, i]() { this->thread_job(i); });

  _finished = false;

  /* start top tree enumeration */
  FT fmaxdisttop = fmaxdist;
  _topenum.enumerate(first, split, fmaxdisttop, fmaxdistexpo, target_coord, subtree, pruning);

  finished = true;
  threadpool.wait_work();
  /* finished enumeration */

  fmaxdist = _bottom_fmaxdist[0];
  for (unsigned i = 0; i < get_threads(); ++i)
    if (_bottom_fmaxdist[i] < fmaxdist)
      fmaxdist = _bottom_fmaxdist[i];

  if (dual && !_evaluator.empty())
  {
    for (auto it = _evaluator.begin(), itend = _evaluator.end(); it != itend; ++it)
      reverse_by_swap(it->second, 0, d - 1);
  }
}

template <typename ZT, typename FT> bool ParallelEnumerationDyn<ZT, FT>::do_work(unsigned i)
{
  vector<enumxt> subtree;
  {
    lockguard lock(_mutex);
    if (_toptrees.empty())
      return false;
    subtree = std::move(_toptrees.back());
    _toptrees.pop_back();
  }
  if (_bottom_fmaxdist[i] != _fmaxdist)
  {
    _bottom_fmaxdist[i] = _fmaxdist;
    _bottom_enums[i].enumerate(_first, _last, _bottom_fmaxdist[i], _fmaxdistexpo, _target_coord,
                               subtree, _pruning, _dual);
  }
  else
    _bottom_enums[i].next_subtree_enumerate(_bottom_fmaxdist[i], _fmaxdistexpo, subtree);
  return true;
}

template <typename ZT, typename FT> void ParallelEnumerationDyn<ZT, FT>::thread_job(unsigned i)
{
  while (true)
  {
    if (!do_work(i))
    {
      if (_finished)
        return;
      else
        std::this_thread::yield();
    }
  }
}

FPLLL_END_NAMESPACE
