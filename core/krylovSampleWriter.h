#pragma once

#include <complex>
#include <cstddef>

namespace TE {

    /**
    * Sink for the sampled wavefunction.
    *
    * The full history of states is Hsize complex numbers per sampling point and
    * cannot be held in memory at the sizes this library targets. The time
    * evolver therefore hands every sample straight to a writer and retains
    * nothing; only the final state is returned in krylovReturn.
    *
    * The implementation lives outside of core, so that core does not depend on
    * HDF5.
    */
    class krylovSampleWriter
    {
    public:
        virtual ~krylovSampleWriter() = default;

        /**
        * Announce the shape of what is about to be written. Called once, before
        * the first sample.
        * @param dim Number of components of a single state
        * @param expectedSamples Number of samples the time evolution intends to take. It is a hint for sizing, not a promise: an evolution that is stopped early writes fewer.
        */
        virtual void beginStates(size_t dim, size_t expectedSamples) = 0;

        /**
        * Consume one sampled state. Called once per sampling point, in order of
        * increasing time.
        * @param state The state at the current sampling point
        * @param dim Number of components, equal to what beginStates announced
        */
        virtual void appendState(const std::complex<double>* state, size_t dim) = 0;

        /**
        * Called once after the last sample. Implementations have to tolerate
        * being called repeatedly and being called without a preceding
        * beginStates, because it is also invoked while an error is on its way
        * out of the time evolution.
        */
        virtual void finishStates() = 0;
    };

}
