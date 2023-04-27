//
// Created by tommaso on 29/03/19.
// Modified by Abhinav and Pietro

#ifndef OPENFPM_PDATA_SUPPORT_HPP
#define OPENFPM_PDATA_SUPPORT_HPP

#include <Space/Shape/Point.hpp>
#include <Vector/vector_dist.hpp>

class Support
{
    // This class is basically a wrapper around a point and a set of offsets.
    // Offsets are stored as they are the data that is mostly required in DCPSE, so we
    // pre-compute and store them, while the positions can be re-computed everytime is
    // necessary (it should almost never be the case) (todo: check if this is required)

private:

    size_t referencePointKey;
    openfpm::vector_std<size_t> keys;

public:

    Support() {};

    Support(const size_t &referencePoint, const openfpm::vector_std<size_t> &keys)
            :referencePointKey(referencePoint),
              keys(keys)
              {}

    Support(const size_t &referencePoint, const std::vector<size_t> &keys)
            :referencePointKey(referencePoint),
              keys(keys.begin(), keys.end())
              {}

    Support(const Support &other)
    : referencePointKey(other.referencePointKey),
      keys(other.keys)
     {}

    size_t size()
    {
        return keys.size();
    }

    const size_t getReferencePointKey() const
    {
        return referencePointKey;
    }

    const openfpm::vector_std<size_t> &getKeys() const
	{
	    return keys;
	}

    openfpm::vector_std<size_t> &getKeys()
	{
	    return keys;
	}

    static bool pack()
    {
        return true;
    }

    static bool packRequest()
    {
        return true;
    }

    template<int ... prp> inline void packRequest(size_t & req) const
    {
        req += sizeof(size_t);
        keys.packRequest(req);
    }

    template<int ... prp> inline void pack(ExtPreAlloc<HeapMemory> & mem, Pack_stat & sts) const
    {
        Packer<size_t,HeapMemory>::pack(mem,referencePointKey,sts);
        keys.template pack<prp ...>(mem,sts);
    }

    template<unsigned int ... prp, typename MemType> inline void unpack(ExtPreAlloc<MemType> & mem, Unpack_stat & ps)
    {
        Unpacker<size_t,MemType>::unpack(mem,referencePointKey,ps);
        keys.template unpack<prp ...>(mem,ps);
    }

};


#endif //OPENFPM_PDATA_SUPPORT_HPP
