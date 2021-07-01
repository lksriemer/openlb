/***************************************************************************
 *   Copyright (C) 2007 by F. P. Beekhof                                   *
 *   fpbeekhof@gmail.com                                                   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <cassert>
#include <limits>
#include <numeric>
#include <cmath>

#include <cvmlcpp/math/Math>
#include <omptl/omptl_algorithm>

namespace cvmlcpp
{


// Adaptor for Rotator3D, for Facet-Normals
// Include dirty hack to avoid const-ness of pairs from maps.
template <typename facet_type,
	typename normal_type = typename facet_type::normal_type,
	typename register_type = typename normal_type::value_type>
class GeoFacetNormRotator
{
	public:
		GeoFacetNormRotator(std::size_t axis, register_type angle) :
			_rotator(axis, angle) { }

		void operator()(const facet_type &f)
		{
			// Write to facet anyway, normal can safely be altered.
			normal_type *p = const_cast<normal_type *>(&f.normal());
			*p = _rotator(f.normal());
// 			f.normal() = _rotator(f.normal());
		}

	private:
		Rotator3D<normal_type> _rotator;
};

// Adaptor for Rotator3D for pair<key, vector>
template <typename Vector, typename register_type = typename Vector::value_type>
class GeoRotator
{
	public:
		GeoRotator(std::size_t axis, register_type angle) :
			_rotator(axis, angle) { }

		const Vector operator()(const Vector &v)
		{ return _rotator(v); }

	private:
		Rotator3D<Vector> _rotator;
};

template <typename T>
Geometry<T>::Geometry()
{
	this->clear();
}

template <typename T>
Geometry<T>::Geometry(const Geometry &that)
{
	*this = that;
}

template <typename T>
Geometry<T> &Geometry<T>::operator=(const Geometry &that)
{
	this->_ptKeygen = that._ptKeygen;
	this->_ftKeygen = that._ftKeygen;

	if (that._minMaxDirty)
		that.recomputeMinMax();
	for (unsigned i = 0u; i < 3u; ++i)
	{
		this->_min[i] = that.min(i);
		this->_max[i] = that.max(i);
	}

	_points.clear();
	_facets.clear();
	_pointKeyMap.clear();
	_facetKeyMap.clear();
	_dirtyFacetNormals.clear();
	_dirtyPointNormals.clear();

	this->_pointNormals  = that._pointNormals;
	this->_pointFacetMap = that._pointFacetMap;

/*
	std::copy(that._pointNormals.begin(), that._pointNormals.end(),
		  this->_pointNormals.begin());
	std::copy(that._pointFacetMap.begin(), that._pointFacetMap.end(),
		  this->_pointFacetMap.begin());
*/
	if(that._facetNormalsDirty)
	{
		assert(that._pointNormalsDirty);
		that.recomputeFacetNormals();
	}
	if (that._pointNormalsDirty)
		that.recomputePointNormals();

	// Copy points
	for (const_pointkey_iterator pkIt = that._pointKeyMap.begin();
	     pkIt != that._pointKeyMap.end(); ++pkIt)
	{
		// Insert point in local set of points
		std::pair<point_iterator, bool> insertResult =
			_points.insert(*pkIt->second);
		assert(insertResult.second);
		// Store local iterator under that same key as in "that"
		_pointKeyMap[pkIt->first] = insertResult.first;
	}

	// Copy facets
	for (const_facetkey_iterator fkIt = that._facetKeyMap.begin();
	     fkIt != that._facetKeyMap.end(); ++fkIt)
	{
		// Insert facet in local set of facets
		std::pair<facet_iterator, bool> insertResult =
			_facets.insert(*fkIt->second);
		assert(insertResult.second);
		// Store local iterator under that same key as in "that"
		_facetKeyMap[fkIt->first] = insertResult.first;
	}

	_facetNormalsDirty = that._facetNormalsDirty;
	_pointNormalsDirty = that._pointNormalsDirty;
	_minMaxDirty	   = that._minMaxDirty;

	return *this;
}

template <typename T>
template <typename SortedPointsIterator, typename PointIDIterator>
Geometry<T>::Geometry(	const SortedPointsIterator ptBegin,
			const SortedPointsIterator ptEnd,
			const PointIDIterator idBegin)
{
	this->loadPoints(ptBegin, ptEnd, idBegin);
}

template <typename T>
template <typename SortedPointsIterator, typename PointIDIterator>
void Geometry<T>::loadPoints(const SortedPointsIterator pointBegin,
			     const SortedPointsIterator pointEnd,
			     const PointIDIterator idBegin)
{
	this->clear();
	assert(_ptKeygen.count() == 0u);

	/*
	 * Create balanced tree of points
	 */
	typedef std::pair<SortedPointsIterator,
			  SortedPointsIterator> PtRange;
	std::vector<PtRange> ptRange[2];
	unsigned e = 0;
	unsigned o = 1;
	ptRange[e].push_back(PtRange(pointBegin, pointEnd));
	while (!ptRange[e].empty())
	{
		assert(ptRange[o].empty());
		for (typename std::vector<PtRange>::const_iterator
		     it = ptRange[e].begin(); it != ptRange[e].end(); ++it)
		{
			const SortedPointsIterator begin1 = it->first;
			const SortedPointsIterator end2   = it->second;

			const std::size_t dist = std::distance(begin1, end2);

			const SortedPointsIterator middle = begin1 + dist / 2u;
			_points.insert(*middle);

			const SortedPointsIterator end1   = middle;
			const SortedPointsIterator begin2 = middle + 1;

			if (begin1 != end1)
				ptRange[o].push_back(PtRange(begin1, end1));
			if (begin2 != end2)
				ptRange[o].push_back(PtRange(begin2, end2));
		}
		ptRange[e].clear();
		using std::swap;
		swap(o, e);
	}

	/*
	 *
	 */
	const std::size_t nPoints = _points.size();
	typedef typename std::iterator_traits<SortedPointsIterator>::difference_type difference_t;
	assert(difference_t(nPoints) == std::distance(pointBegin, pointEnd));

	std::vector<point_iterator> pointIterators(nPoints);
	std::size_t index = 0u;
	for (point_iterator pIt = _points.begin(); pIt != _points.end(); ++pIt)
		pointIterators[index++] = pIt;


	// Generate and sort ID's for points. Points and ID's are now aligned
	PointIDIterator idEnd = idBegin;
	for (std::size_t i = 0; i < nPoints; ++i, ++idEnd)
		*idEnd = _ptKeygen();
	omptl::sort(idBegin, idEnd);

	/*
	 * Create balanced tree of iterators-to-points
	 */
	typedef std::pair<PointIDIterator, PointIDIterator> IDRange;
	std::vector<IDRange> idRange[2];
	e = 0;
	o = 1;
	idRange[e].push_back(IDRange(idBegin, idEnd));
	while (!idRange[e].empty())
	{
		assert(idRange[o].empty());
		for (typename std::vector<IDRange>::const_iterator
		     it = idRange[e].begin(); it != idRange[e].end(); ++it)
		{
			const PointIDIterator begin1 = it->first;
			const PointIDIterator end2   = it->second;

			const std::size_t dist   = std::distance(begin1, end2);
			const PointIDIterator middle = begin1 + dist / 2u;

			const std::size_t offset = std::distance(idBegin, middle);
			_pointKeyMap[*middle] = pointIterators[offset];
			_dirtyPointNormals.insert(*middle);

			const PointIDIterator end1   = middle;
			const PointIDIterator begin2 = middle + 1;

			if (begin1 != end1)
				idRange[o].push_back(IDRange(begin1, end1));
			if (begin2 != end2)
				idRange[o].push_back(IDRange(begin2, end2));
		}
		idRange[e].clear();
		using std::swap;
		swap(o, e);
	}

	_minMaxDirty = true;
	this->setNormalsDirty();

	assert(_ptKeygen.count() == this->nrPoints());
}

template <typename T>
bool Geometry<T>::operator==(const Geometry &that) const
{
	if ( (this->nrFacets() != that.nrFacets()) ||
	     (this->nrPoints() != that.nrPoints()) )
		return false;

	for (std::size_t i = 0u; i < 3u; ++i)
		if ( (this->min(i) != that.min(i)) ||
		     (this->max(i) != that.max(i)) )
			return false;

	// Compare points.
	assert(this->_pointKeyMap.size() == that._pointKeyMap.size());
	for (const_pointkey_iterator pkItThis = this->_pointKeyMap.begin(),
	     pkItThat = that. _pointKeyMap.begin();
	     pkItThat != that._pointKeyMap.end(); ++pkItThis, ++pkItThat)
		if ( (pkItThis->first != pkItThat->first) ||
		     (*pkItThis->second != *pkItThat->second) )
			return false;

	assert(this->_facetKeyMap.size() == that._facetKeyMap.size());
	for (const_facetkey_iterator fkItThis = this->_facetKeyMap.begin(),
	     fkItThat = that. _facetKeyMap.begin();
	     fkItThat != that._facetKeyMap.end(); ++fkItThis, ++fkItThat)
		if ( ( fkItThis->first  !=  fkItThat->first) ||
		     (*fkItThis->second != *fkItThat->second) )
			return false;

	return true;
}

template <typename T>
Geometry<T>::~Geometry()
{
// 	this->clear();
}

template <typename T>
void Geometry<T>::clear()
{
	_ptKeygen.reset();
	_ftKeygen.reset();

	std::fill(_min.begin(), _min.end(),  std::numeric_limits<T>::max());
	std::fill(_max.begin(), _max.end(), -std::numeric_limits<T>::max());

	_points.clear();
	_facets.clear();
	_pointKeyMap.clear();
	_facetKeyMap.clear();

	_pointNormals.clear();
	_pointFacetMap.clear();

	_dirtyFacetNormals.clear();
	_dirtyPointNormals.clear();

	_facetNormalsDirty = false;
	_pointNormalsDirty = false;
	_minMaxDirty	   = false;
}

template <typename T>
const typename Geometry<T>::facet_type &Geometry<T>::
facet(const std::size_t key) const
{
	// Internal Consistency check
	assert(_facetKeyMap.size() == _facets.size());

	// User has asked for non-existent key ?
	assert(_facetKeyMap.find(key) != _facetKeyMap.end());

	return *_facetKeyMap[key];
}

template <typename T>
const typename Geometry<T>::vector_type &Geometry<T>::
facetNormal(const std::size_t key) const
{
	if(_facetNormalsDirty)
	{
		assert(_pointNormalsDirty);
		this->recomputeFacetNormals();
	}

	return this->facet(key).normal();
}

template <typename T>
const typename Geometry<T>::vector_type &Geometry<T>::
pointNormal(const std::size_t key) const
{
	// Internal Consistency check
	assert(_pointNormals.size() == _points.size());

	// User has asked for non-existent key ?
	assert(_pointNormals.find(key) != _pointNormals.end());

	if (_pointNormalsDirty)
		this->recomputePointNormals();

	return _pointNormals[key];
}

template <typename T>
const typename Geometry<T>::point_type &Geometry<T>::
point(const std::size_t key) const
{
	// Internal Consistency check
	assert(_pointKeyMap.size() == _points.size());

	// User has asked for non-existent key ?
	assert(_pointKeyMap.find(key) != _pointKeyMap.end());
	return *_pointKeyMap[key];
}

template <typename T>
const std::set<std::size_t> &Geometry<T>::
facetsHavingPoint(const std::size_t key) const
{
	assert(_pointFacetMap.find(key) != _pointFacetMap.end());
	return _pointFacetMap[key];
}

template <typename T>
std::size_t Geometry<T>::addPoint(const value_type &x,
				  const value_type &y,
				  const value_type &z)
{
	return this->addPoint(point_type(x, y, z));
}

template <typename T>
std::size_t Geometry<T>::addPoint(const point_type &point)
{
	assert(point[X] >= value_type(0.0) || point[X] <= value_type(0.0));
	assert(point[Y] >= value_type(0.0) || point[Y] <= value_type(0.0));
	assert(point[Z] >= value_type(0.0) || point[Z] <= value_type(0.0));

	std::size_t key;
	if (_points.empty())
	{
		// It is really a new point.
		key = _ptKeygen.generate();
		assert(_pointKeyMap.find(key) == _pointKeyMap.end());
		_pointKeyMap[key] = _points.insert(point).first;
		std::fill(_pointNormals[key].begin(),
			  _pointNormals[key].end(), value_type(0));
	}
	else
	{
		const const_point_iterator pIt = _points.lower_bound(point);

		// Does this point exist already ?
		if ( (pIt != _points.end()) && (*pIt == point) )
		{
			for (const_pointkey_iterator pkIt= _pointKeyMap.begin();
			     pkIt != _pointKeyMap.end(); ++pkIt)
				if (pkIt->second == pIt)
					return pkIt->first;
			assert(false);
		}

		// It is really a new point.
		key = _ptKeygen.generate();
		assert(_pointKeyMap.find(key) == _pointKeyMap.end());
		_pointKeyMap[key] = _points.insert(pIt, point);
		std::fill(_pointNormals[key].begin(),
			_pointNormals[key].end(), value_type(0));
	}

	// If dirty, there's no point in updating
	if (!_minMaxDirty)
	{
		for (unsigned i = 0; i < 3u; ++i)
			_min[i] = std::min(_min[i], point[i]);
		for (unsigned i = 0; i < 3u; ++i)
			_max[i] = std::max(_max[i], point[i]);
	}

	return key;
}

template <typename T>
bool Geometry<T>::updatePoint(const std::size_t key, const value_type &x,
			      const value_type &y,   const value_type &z)
{
	return this->updatePoint(key, point_type(x, y, z));
}

template <typename T>
bool Geometry<T>::updatePoint(const std::size_t key, const point_type &p)
{
	// Does this point really exist ?
	pointkey_iterator pkIt = _pointKeyMap.find(key);
	if (pkIt == _pointKeyMap.end())
		return false;

	// Erase old point
	_points.erase(pkIt->second);
	_minMaxDirty = true;

	// Insert new point, update keys
	std::pair<point_iterator, bool> insertResult = _points.insert(p);
	if (!insertResult.second)
		return false;
	pkIt->second = insertResult.first;

	// Mark involved normals to be recomputed
	_dirtyPointNormals.insert(key);
	_dirtyFacetNormals.insert(_pointFacetMap[key].begin(),
				  _pointFacetMap[key].end());

	this->setNormalsDirty();

	return true;
}

template <typename T>
bool Geometry<T>::erasePoint(const std::size_t key)
{
	const pointkey_iterator pkIt = _pointKeyMap.find(key);
	if (pkIt == _pointKeyMap.end())
		return false;
	if (_pointFacetMap[key].size() != 0u)
		return false;

	_points.erase(pkIt->second);
	_pointKeyMap.erase(pkIt);
	_pointFacetMap.erase(key);
	_pointNormals.erase(key);
	_dirtyPointNormals.erase(key);

	_minMaxDirty = true;
	this->setNormalsDirty();

	return true;
}

template <typename T>
typename Geometry<T>::const_facet_iterator Geometry<T>::facetsBegin() const
{
	if(_facetNormalsDirty)
		this->recomputeFacetNormals();
	return _facets.begin();
}

template <typename T>
typename Geometry<T>::const_facet_iterator Geometry<T>::facetsEnd() const
{
	return _facets.end();
}

template <typename T>
std::size_t Geometry<T>::addFacet(const std::size_t a,
				const std::size_t b, const std::size_t c)
{
	const std::size_t fpts [] = {a, b, c};
	return this->addFacet(facet_type(fpts, fpts+3));
}

template <typename T>
std::size_t Geometry<T>::addFacet(const facet_type &facet)
{
	// User input ok ? I.e., do the points exist ?
	assert(_pointKeyMap.find(facet[A]) != _pointKeyMap.end());
	assert(_pointKeyMap.find(facet[B]) != _pointKeyMap.end());
	assert(_pointKeyMap.find(facet[C]) != _pointKeyMap.end());

	std::size_t key;

	if (_facets.empty())
	{
		key = _ftKeygen();
		assert(_facetKeyMap.find(key) == _facetKeyMap.end());
		_facetKeyMap[key] = _facets.insert(facet).first;
	}
	else
	{
		const facet_iterator fIt = _facets.lower_bound(facet);

		// Does this facet exist already ?
		if ( (fIt != _facets.end()) && (*fIt == facet) )
		{
			for (const_facetkey_iterator fkIt= _facetKeyMap.begin();
			     fkIt != _facetKeyMap.end(); ++fkIt)
				if (fkIt->second == fIt)
					return fkIt->first;
			assert(false);
		}

		// It is really a new facet.
		key = _ftKeygen();
		assert(_facetKeyMap.find(key) == _facetKeyMap.end());
		_facetKeyMap[key] = _facets.insert(fIt, facet);
	}

	_dirtyPointNormals.insert(facet.begin(), facet.end());
	_pointNormalsDirty = true;

	if (facet.normal() == T(0))
	{
		_dirtyFacetNormals.insert(key);
		_facetNormalsDirty = true;
	}

	_pointFacetMap[facet[A]].insert(key);
	_pointFacetMap[facet[B]].insert(key);
	_pointFacetMap[facet[C]].insert(key);

	return key;
}

template <typename T>
template <typename SortedFacetsIterator,
		typename FacetIDIterator>
void Geometry<T>::loadFacets(const SortedFacetsIterator facetBegin,
			     const SortedFacetsIterator facetEnd,
			     const FacetIDIterator idBegin)
{
	// Clean out old facet data
	_ftKeygen.reset();
	_facets.clear();
	_facetKeyMap.clear();
	_pointFacetMap.clear();
	this->setNormalsDirty();

	/*
	 * Create balanced tree of facets
	 */
	typedef std::pair<SortedFacetsIterator,
			  SortedFacetsIterator> FRange;
	std::vector<FRange> fRange[2];
	unsigned e = 0;
	unsigned o = 1 - e;
	fRange[e].push_back(FRange(facetBegin, facetEnd));
	while (!fRange[e].empty())
	{
		assert(fRange[o].empty());
		for (typename std::vector<FRange>::const_iterator
		     it = fRange[e].begin(); it != fRange[e].end(); ++it)
		{
			const SortedFacetsIterator begin1 = it->first;
			const SortedFacetsIterator end2   = it->second;

			const std::size_t dist = std::distance(begin1, end2);

			const SortedFacetsIterator middle = begin1 + dist / 2u;
			_facets.insert(*middle);
			if (_dirtyPointNormals.size() < this->nrPoints())
			{
				_dirtyPointNormals.insert((*middle)[A]);
				_dirtyPointNormals.insert((*middle)[B]);
				_dirtyPointNormals.insert((*middle)[C]);
			}

			const SortedFacetsIterator end1   = middle;
			const SortedFacetsIterator begin2 = middle + 1;

			if (begin1 != end1)
				fRange[o].push_back(FRange(begin1, end1));
			if (begin2 != end2)
				fRange[o].push_back(FRange(begin2, end2));
		}
		fRange[e].clear();
		using std::swap;
		swap(o, e);
	}

	/*
	 *
	 */
	const std::size_t nFacets = _facets.size();
	std::vector<facet_iterator> facetIterators(nFacets);
	std::size_t index = 0u;
	for (facet_iterator fIt = _facets.begin(); fIt != _facets.end(); ++fIt)
	{
		facetIterators[index++] = fIt;
		_dirtyPointNormals.insert(fIt->begin(), fIt->end());
	}


	// Generate and sort ID's for points. Points and ID's are now aligned
	FacetIDIterator idEnd = idBegin;
	for (std::size_t i = 0; i < nFacets; ++i, ++idEnd)
		*idEnd = _ftKeygen();
	omptl::sort(idBegin, idEnd);

	/*
	 * Create balanced tree of iterators-to-facets
	 */
	typedef std::pair<FacetIDIterator, FacetIDIterator> IDRange;
	std::vector<IDRange> idRange[2];
	e = 0;
	o = 1;
	idRange[e].push_back(IDRange(idBegin, idEnd));
	while (!idRange[e].empty())
	{
		assert(idRange[o].empty());
		for (typename std::vector<IDRange>::const_iterator
		     it = idRange[e].begin(); it != idRange[e].end(); ++it)
		{
			const FacetIDIterator begin1 = it->first;
			const FacetIDIterator end2   = it->second;

			const std::size_t dist   = std::distance(begin1, end2);
			const FacetIDIterator middle = begin1 + dist / 2u;

			const std::size_t offset = std::distance(idBegin, middle);
			_facetKeyMap[*middle] = facetIterators[offset];

			for (std::size_t p = 0; p < 3; ++p)
				_pointFacetMap[(*facetIterators[offset])[p]]
						.insert(*middle);

			const FacetIDIterator end1   = middle;
			const FacetIDIterator begin2 = middle + 1;

			if (begin1 != end1)
				idRange[o].push_back(IDRange(begin1, end1));
			if (begin2 != end2)
				idRange[o].push_back(IDRange(begin2, end2));
		}
		idRange[e].clear();
		using std::swap;
		swap(o, e);
	}

	assert(_ftKeygen.count() == this->nrFacets());
	assert(nFacets == this->nrFacets());
}

template <typename T>
bool Geometry<T>::updateFacet(	const std::size_t key, const std::size_t a,
				const std::size_t b, const std::size_t c)
{
	const std::size_t fpts [] = {a, b, c};
	const facet_type facet = facet_type(fpts, fpts+3);
	return this->updateFacet(key, facet);
}

template <typename T>
bool Geometry<T>::updateFacet(const std::size_t key, const facet_type &facet)
{
	// Does this facet really exist ?
	facetkey_iterator fkIt = _facetKeyMap.find(key);
	if (fkIt == _facetKeyMap.end())
		return false;

	// Do the points in the new facet exist ?
	if ( (_pointKeyMap.find(facet[A]) == _pointKeyMap.end()) ||
	     (_pointKeyMap.find(facet[B]) == _pointKeyMap.end()) ||
	     (_pointKeyMap.find(facet[C]) == _pointKeyMap.end()) )
		return false;

	// Erase old facet
	_pointFacetMap[(*fkIt->second)[A]].erase(key);
	_pointFacetMap[(*fkIt->second)[B]].erase(key);
	_pointFacetMap[(*fkIt->second)[C]].erase(key);
	_facets.erase(fkIt->second);

	// Insert new facet, update keys
	std::pair<facet_iterator, bool> insertResult = _facets.insert(facet);
	if (!insertResult.second)
		return false;
	fkIt->second = insertResult.first;

	// Mark involved normals to be recomputed
	_dirtyPointNormals.insert(fkIt->second->begin(), fkIt->second->end());
	_dirtyFacetNormals.insert(key);

	_minMaxDirty = true;
	this->setNormalsDirty();

	return true;
}

template <typename T>
bool Geometry<T>::eraseFacet(const std::size_t key)
{
	const facetkey_iterator fkIt = _facetKeyMap.find(key);
	if (fkIt == _facetKeyMap.end())
		return false;

	for (typename facet_type::const_triangle_iterator fIt =
	     fkIt->second->begin(); fIt != fkIt->second->end(); ++fIt)
		_pointFacetMap[*fIt].erase(key);

	_dirtyFacetNormals.erase(key);
	_facets.erase(fkIt->second);
	_facetKeyMap.erase(fkIt);

	this->setNormalsDirty();


	return true;
}

template <typename T>
T Geometry<T>::min(const unsigned dim) const
{
	// Member function should not be called without point
	assert(this->nrPoints() > 0);

	if (_minMaxDirty)
		this->recomputeMinMax();

	return _min[dim];
}

template <typename T>
T Geometry<T>::max(const unsigned dim) const
{
	// Member function should not be called without point
	assert(this->nrPoints() > 0);

	if (_minMaxDirty)
		this->recomputeMinMax();

	return _max[dim];
}

template <typename T>
void Geometry<T>::center()
{
	const T dx = (this->max(X) + this->min(X)) * 0.5;
	const T dy = (this->max(Y) + this->min(Y)) * 0.5;
	const T dz = (this->max(Z) + this->min(Z)) * 0.5;

	this->translate(-dx, -dy, -dz);
}

template <typename T>
void Geometry<T>::centerMass()
{
	typename Geometry<T>::point_type mass =
		std::accumulate(this->pointsBegin(),
		 	        this->pointsEnd(), point_type(0.0));

	this->translate(mass * T(-1));
}

template <typename T>
void Geometry<T>::scale(const T factor)
{
	if (this->nrPoints() == 0)
		return;

	// Bogus user input ?
	assert(factor > 0.0);

	// Very dirty, but should be safe, provided factor > 0
	for (const_point_iterator it = this->pointsBegin();
	     it != this->pointsEnd(); ++it)
	{
		point_type * p = const_cast<point_type *>(&(*it));
		*p = *it * factor;
	}

/*
	std::map<std::size_t, point_type> tempPts;

	// Scale points and store in temporary map
	for (const_pointkey_iterator it = _pointKeyMap.begin();
	     it != _pointKeyMap.end(); ++it)
		tempPts[it->first] = *(it->second) * factor;

	// Delete old set of points
	_points.clear();

	// Insert new points and update keymap
	typedef typename std::map<std::size_t, point_type>::const_iterator TmpIt;
	for (TmpIt it = tempPts.begin(); it != tempPts.end(); ++it)
	{
		// Insert
		std::pair<const_point_iterator, bool>
			insResult = _points.insert(it->second);

		// Should be impossible
		assert(insResult.second);

		// Update keymap
		_pointKeyMap[it->first] = insResult.first;
	}
*/
	// Fix min & max
	if (!_minMaxDirty)
	{
		_min *= factor;
		_max *= factor;
	}
}

template <typename T>
void Geometry<T>::scaleTo(const T maxLength)
{
	const T dx = this->max(X) - this->min(X);
	const T dy = this->max(Y) - this->min(Y);
	const T dz = this->max(Z) - this->min(Z);

	const T maxLen = std::max( dx, std::max(dy, dz) );
	this->scale(maxLength / maxLen);
}

template <typename T>
template <class Vector>
void Geometry<T>::translate(const Vector &v)
{
	if (this->nrPoints() == 0)
		return;

	// Works only for floating-point types!
	const T e = T(2.0) * std::numeric_limits<T>::min();
	if (!( (std::abs(v[X])>e) || (std::abs(v[Y])>e) || (std::abs(v[Z])>e) ))
		return;

	typedef std::pair<std::size_t, point_type> pPair;
	std::vector<pPair> tempPts;

	// Translate points and store in temporary map
	for (const_pointkey_iterator it = _pointKeyMap.begin();
	     it != _pointKeyMap.end(); ++it)
		tempPts.push_back(pPair(it->first, *it->second + v));

	// Delete old set of points
	assert(_ptKeygen.count() >= this->nrPoints());
	_points.clear();
	_pointKeyMap.clear();

	// Insert new points and update keymap
	typedef typename std::vector<pPair>::const_iterator TmpPIt;
	typedef PairFirstCompare<std::size_t, point_type,
			std::less<std::size_t> > PairLess;

	KeyGenerator<std::size_t> reGen;
	for (std::size_t i = 0u; i < _ptKeygen.count(); ++i)
	{
		const std::size_t key = reGen.generate();
		const TmpPIt it = std::lower_bound(tempPts.begin(),
						tempPts.end(), key, PairLess());

		if ( (it != tempPts.end()) && (it->first == key) )
		{
			// Insert
			std::pair<const_point_iterator, bool>
				insResult = _points.insert(it->second);
			assert(insResult.second); // Should be impossible

			// Update keymap
			_pointKeyMap.insert(std::pair<std::size_t, point_iterator>
						(it->first, insResult.first));
		}
	}
	assert(reGen.count() == _ptKeygen.count());

	// Fix min & max
	if (!_minMaxDirty)
	{
		_min += v;
		_max += v;
	}
}

template <typename T>
void Geometry<T>::translate(const T dx, const T dy, const T dz)
{
	this->translate(vector_type(dx, dy, dz));
}

template <typename T>
void Geometry<T>::rotate(const std::size_t axis, const T angle)
{
	// Used supplied invalid axis ?
	assert( (axis == X) || (axis == Y) || (axis == Z) );

	if (this->nrPoints() == 0)
		return;

	// Works only for floating-point types!
	const T e = T(2.0) * std::numeric_limits<T>::min();
	if ( std::abs(std::fmod(angle, T(2.0)*Constants<T>::pi())) < e )
		return;

	/*
	 * Rotate Points
	 */

	// Rotate and place in temporary map
	typedef std::pair<std::size_t, point_type> pPair;
	std::vector<pPair> tempPts;

	GeoRotator<point_type> grp((axis), angle);

	for (const_pointkey_iterator it = _pointKeyMap.begin();
	     it != _pointKeyMap.end(); ++it)
		tempPts.push_back(pPair(it->first, grp(*(it->second))));

	// Insert new points and update keymap
	_points.clear();
	_pointKeyMap.clear();
	KeyGenerator<std::size_t> reGen;
	for (std::size_t i = 0u; i < _ptKeygen.count(); ++i)
	{
		typedef typename std::vector<pPair>::const_iterator TmpPIt;
		typedef PairFirstCompare<std::size_t, point_type,
						std::less<std::size_t> > PairLess;

		const std::size_t key = reGen.generate();
		const TmpPIt it = std::lower_bound(tempPts.begin(),
						tempPts.end(), key, PairLess());
		if ( (it != tempPts.end()) && (it->first == key) )
		{
			// Insert
			std::pair<const_point_iterator, bool>
				insResult = _points.insert(it->second);
			assert(insResult.second); // Should be impossible

			// Update keymap
			_pointKeyMap.insert(std::pair<std::size_t, point_iterator>
						(it->first, insResult.first));
		}
	}

	/*
	 * Rotate Normals
	 */
	std::for_each(_facets.begin(), _facets.end(),
			GeoFacetNormRotator<facet_type>(axis, angle));

	/*
	 * Rotate Point-Normals
	 */
	std::map<std::size_t, vector_type> newPtNormals;
	typedef MapPairOperateInserter<
				std::map<std::size_t, vector_type>,
				GeoRotator<vector_type> > PointNormsInserter;


	GeoRotator<vector_type> grv(axis, angle);
	PointNormsInserter pni((newPtNormals), grv);
	std::for_each(_pointNormals.begin(), _pointNormals.end(), pni);

	_pointNormals.swap(newPtNormals);

	// Didn't fix min/max
	_minMaxDirty = true;
}

template <typename T>
void Geometry<T>::recomputeMinMax() const
{
	std::fill(_min.begin(), _min.end(),  std::numeric_limits<T>::max());
	std::fill(_max.begin(), _max.end(), -std::numeric_limits<T>::max());

	for (typename std::set<point_type>::const_iterator p = _points.begin();
	     p != _points.end(); ++p)
	{
		for (std::size_t i = 0; i < 3u; ++i)
			_min[i] = std::min(_min[i], (*p)[i]);
		for (std::size_t i = 0; i < 3u; ++i)
			_max[i] = std::max(_max[i], (*p)[i]);
	}

	_minMaxDirty = false;
}

template <typename T>
void Geometry<T>::recomputeFacetNormals() const
{
	// Recompute facet normals
	for (std::set<std::size_t>::const_iterator
	     it = _dirtyFacetNormals.begin();
	     it != _dirtyFacetNormals.end(); ++it)
	{
// std::cout << "recomputeFacetNormals() Facet: " << *it << std::endl;
		facet_type fac = *_facetKeyMap[*it];

		// Valid Facet ?
		assert(fac[A] != fac[B]);
		assert(fac[A] != fac[C]);
		assert(fac[B] != fac[C]);
// std::cout << "recomputeFacetNormals() abc: " << fac[A] << " " << fac[B]
// 	<< " " << fac[C] << std::endl;
		// Points exist ?
		assert(_pointKeyMap.find(fac[A]) != _pointKeyMap.end());
		assert(_pointKeyMap.find(fac[B]) != _pointKeyMap.end());
		assert(_pointKeyMap.find(fac[C]) != _pointKeyMap.end());

		// Changing the facet normal will change the point normals
		_dirtyPointNormals.insert(fac[A]);
		_dirtyPointNormals.insert(fac[B]);
		_dirtyPointNormals.insert(fac[C]);

		const vector_type ab =
				*_pointKeyMap[fac[B]] - *_pointKeyMap[fac[A]];
		const vector_type ac =
				*_pointKeyMap[fac[C]] - *_pointKeyMap[fac[A]];
		assert(dotProduct(ab, ab) > 0.0f);
		assert(dotProduct(ac, ac) > 0.0f);

		const vector_type norm = crossProduct(ac, ab);
		if (!(dotProduct(norm, norm) > 0.0f))
		{
			// Computation of normal failed! Bugger. Now what ?
			facet_type *p = const_cast<facet_type *>
						(&(*_facetKeyMap[*it]));
			p->normal() = vector_type(value_type(0.0));
			continue;
		}
		assert(dotProduct(norm, norm) > 0.0f);
		facet_type *p = const_cast<facet_type *>(&(*_facetKeyMap[*it]));
		p->normal() = norm / std::sqrt(dotProduct(norm, norm));
	}

	_dirtyFacetNormals.clear();
	_facetNormalsDirty = false;
}

template <typename T>
void Geometry<T>::recomputePointNormals() const
{
	if(_facetNormalsDirty)
		this->recomputeFacetNormals();

	// Recompute point normals
	for (std::set<std::size_t>::const_iterator
	     it = _dirtyPointNormals.begin();
	     it != _dirtyPointNormals.end(); ++it)
	{
		vector_type pNormal(value_type(0.0));

		for (std::set<std::size_t>::const_iterator
		     normIt = _pointFacetMap[*it].begin();
		     normIt != _pointFacetMap[*it].end(); ++normIt)
			pNormal += _facetKeyMap[*normIt]->normal();
// 		pNormal /= value_type(_pointFacetMap[*it].size());
		pNormal /= std::sqrt(dotProduct(pNormal, pNormal));

		_pointNormals[*it] = pNormal;
	}

	_dirtyPointNormals.clear();
	_pointNormalsDirty = false;
}

template <typename T>
inline void Geometry<T>::setNormalsDirty()
{
	_facetNormalsDirty = true;
	_pointNormalsDirty = true;
}

} // namespace cvmlcpp
