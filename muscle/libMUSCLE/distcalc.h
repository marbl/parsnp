#ifndef DistCalc_h
#define DistCalc_h

namespace muscle {

typedef float dist_t;
const dist_t BIG_DIST = (dist_t) 1e29;

class DistFunc;

class DistCalc
	{
public:
	virtual void CalcDistRange(unsigned i, dist_t Dist[]) const = 0;
	virtual unsigned GetCount() const = 0;
	virtual unsigned GetId(unsigned i) const = 0;
	virtual const char *GetName(unsigned i) const = 0;
	};

class DistCalcDF : public DistCalc
	{
public:
	void Init(const DistFunc &DF);
	virtual void CalcDistRange(unsigned i, dist_t Dist[]) const;
	virtual unsigned GetCount() const;
	virtual unsigned GetId(unsigned i) const;
	virtual const char *GetName(unsigned i) const;

private:
	const DistFunc *m_ptrDF;
	};

class DistCalcMSA : public DistCalc
	{
public:
	void Init(const MSA &msa, DISTANCE Distance);
	virtual void CalcDistRange(unsigned i, dist_t Dist[]) const;
	virtual unsigned GetCount() const;
	virtual unsigned GetId(unsigned i) const;
	virtual const char *GetName(unsigned i) const;

private:
	const MSA *m_ptrMSA;
	DISTANCE m_Distance;
	};

} // namespace muscle

#endif	// DistCalc_h
