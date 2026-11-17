CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17
CPPFLAGS += -I include/
LDLIBS += -lgmp

SRCS = count.cpp rand.cpp counthelp.cpp

# test.dnf has exactly 43675/2^20 = 0.04165172... satisfying assignments
# (computed by inclusion-exclusion). count uses eps=delta=0.05, so each run is
# within relative error eps of this value with probability >= 0.95. The test
# takes the median of 5 independent runs, which fails spuriously only if 3 of
# the 5 runs miss (probability ~0.1%).
EXPECTED = 0.0416517257690430
EPS = 0.05
RUNS = 5

all: count

count: $(SRCS) include/DMCrand.h include/DMChelp.h include/bit_array.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $(SRCS) $(LDLIBS)

test: count
	@for i in $$(seq $(RUNS)); do ./count; done | sort -g | \
	awk -v ref=$(EXPECTED) -v eps=$(EPS) -v runs=$(RUNS) ' \
	    { v[NR] = $$1; printf "run: %s\n", $$1 } \
	    END { \
	        if (NR != runs) { print "FAIL: expected " runs " estimates, got " NR; exit 1 } \
	        med = v[int((NR + 1) / 2)]; \
	        err = (med - ref) / ref; if (err < 0) err = -err; \
	        printf "median %s, exact %s, relative error %.4f (eps %s)\n", med, ref, err, eps; \
	        if (err > eps) { print "FAIL"; exit 1 } \
	        print "PASS" \
	    }'

clean:
	rm -f count

.PHONY: all test clean
