//
// Created by tommaso on 29/03/19.
// Modified by Abhinav and Pietro

#ifndef OPENFPM_PDATA_SUPPORT_HPP
#define OPENFPM_PDATA_SUPPORT_HPP

#include <Space/Shape/Point.hpp>
#include <Vector/vector_dist.hpp>

template<typename vector_type>
class Support
{
    // This class is basically a wrapper around a point and a set of offsets.
    // Offsets are stored as they are the data that is mostly required in DCPSE, so we
    // pre-compute and store them, while the positions can be re-computed everytime is
    // necessary (it should almost never be the case) (todo: check if this is required)

private:
    const vector_type &domain;
    const size_t referencePointKey;
    const std::vector<size_t> keys;
    std::vector<Point<vector_type::dims,typename vector_type::stype>> offsets;

public:
    Support() {};

    Support(const vector_type &domain, const size_t &referencePoint, const std::vector<size_t> &keys)
            : domain(domain),
              referencePointKey(referencePoint),
              keys(keys),
              offsets(computeOffsets(referencePoint, keys)) {}

    Support(const Support<vector_type> &other);

    size_t size();

    const Point<vector_type::dims,typename vector_type::stype> getReferencePoint() const;

    const size_t getReferencePointKey() const;

    const std::vector<size_t> &getKeys() const;

    const std::vector<Point<vector_type::dims,typename vector_type::stype>> &getOffsets() const;

    std::vector<Point<vector_type::dims,typename vector_type::stype>>
    RecomputeOffsets();

private:
    std::vector<Point<vector_type::dims,typename vector_type::stype>>
    computeOffsets(const size_t referencePoint, const std::vector<size_t> &keys);
};

template<typename vector_type>
std::vector<Point<vector_type::dims,typename vector_type::stype>>
Support<vector_type>::computeOffsets(const size_t referencePoint, const std::vector<size_t> &keys)
{
    std::vector<Point<vector_type::dims,typename vector_type::stype>> offsets;
    for (auto &otherK : keys)
    {
        Point<vector_type::dims,typename vector_type::stype> curOffset(domain.getPos(referencePoint));
        curOffset -= domain.getPos(otherK);
        offsets.push_back(curOffset);
    }
    return offsets;
}

template<typename vector_type>
std::vector<Point<vector_type::dims,typename vector_type::stype>>
Support<vector_type>::RecomputeOffsets()
{
    offsets.clear();
    for (auto &otherK : keys)
    {
        Point<vector_type::dims,typename vector_type::stype> curOffset(domain.getPos(referencePointKey));
        curOffset -= domain.getPos(otherK);
        offsets.push_back(curOffset);
    }
    return offsets;
}

template<typename vector_type>
const Point<vector_type::dims,typename vector_type::stype> Support<vector_type>::getReferencePoint() const
{
    return Point<vector_type::dims,typename vector_type::stype>(domain.getPos(referencePointKey));
}

template<typename vector_type>
const std::vector<Point<vector_type::dims,typename vector_type::stype>> &Support<vector_type>::getOffsets() const
{
    return offsets;
}

template<typename vector_type>
size_t Support<vector_type>::size()
{
    return offsets.size();
}

template<typename vector_type>
Support<vector_type>::Support(const Support<vector_type> &other)
        : domain(other.domain),
          referencePointKey(other.referencePointKey),
          keys(other.keys),
          offsets(other.offsets) {}

template<typename vector_type>
const size_t Support<vector_type>::getReferencePointKey() const
{
    return referencePointKey;
}

template<typename vector_type>
const std::vector<size_t> &Support<vector_type>::getKeys() const
{
    return keys;
}

#endif //OPENFPM_PDATA_SUPPORT_HPP
