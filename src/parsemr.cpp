#include "ExternalSequence.h"
#include <stdio.h>

//#include <chrono> // for std::this_thread::sleep_for()
//#include <thread> // for std::this_thread::sleep_for()
//using namespace std::chrono_literals;

// optionally use SEQ_NAMESPACE
#ifdef SEQ_NAMESPACE
using namespace SEQ_NAMESPACE;
#endif

void print_usage()
{
    printf("parsemr: standalone Pulseq file loader\n");
    printf("         This is a small C++ application intended \n"
           "         to test and demo the Pulseq C++ library code\n");
    printf("usage: parsemr <sequence_file>\n");
}

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        print_usage();
        return 1;
    }

    std::string pulseqFilePath = argv[1];
    // TODO: check if file exists and whether this is a Pulseq file

    // we create a dummy scope to see later on that the memoty has been freed without issues
    {
        // create an external sequence object and load the file
        ExternalSequence SequenceData;
        bool             bSeqDataPrepared = SequenceData.load(pulseqFilePath);

        if (!bSeqDataPrepared)
        {
            printf("Failed to load the sequence \'%s\', exiting...", pulseqFilePath.c_str());
            return 1;
        }

        // check the signature
        if (SequenceData.isSigned())
        {
            std::cout << "The file contains " << SequenceData.getSignatureType() << " signature "
                      << SequenceData.getSignature() << std::endl;
            std::cout << "Signature check " << (SequenceData.isSignatureCheckSucceeded() ? "succeeded" : "FAILED")
                      << std::endl;
        }
        else
            std::cout << "The file contains no signature" << std::endl;

        // list all definitions
        std::vector<std::string> vecDefinitionKeys = SequenceData.GetAllDefinitions();
        if (vecDefinitionKeys.empty())
            std::cout << "The file contains no definitions" << std::endl;
        else
        {
            std::cout << "The file contains following definitions:" << std::endl;
            for (int i = 0; i < vecDefinitionKeys.size(); ++i)
            {
                std::cout << vecDefinitionKeys[i] << " : \t";
                std::cout << SequenceData.GetDefinitionStr(vecDefinitionKeys[i]) << std::endl;
            }
        }
        std::cout << std::endl;

        // scan through the sequence, get blocks
        int                      nRF             = 0;
        int                      nGradArb        = 0;
        int                      nGradTrap       = 0;
        int                      nGradExtTrap    = 0;
        int                      nLBL            = 0;
        int                      nADC            = 0;
        int                      nSoDe           = 0;
        int                      nRot            = 0;
        int64_t                  llTotalDuration = 0;
        LabelStateAndBookkeeping labelStateAndBookkeeping;
        int                      iB;
        for (iB = 0; iB < SequenceData.GetNumberOfBlocks(); ++iB)
        {
            SeqBlock* pBlock = SequenceData.GetBlock(iB);
            // look into the block, unpack content and collect some data or statistics

            if (pBlock->isRF())
            {
                ++nRF;
            }

            // Process each gradient channel
            for (int iChannel = 0; iChannel < NUM_GRADS; iChannel++)
            {
                if (pBlock->isArbitraryGradient(iChannel))
                {
                    ++nGradArb;
                }
                else if (pBlock->isExtTrapGradient(iChannel))
                {
                    ++nGradExtTrap;
                }
                else if (pBlock->isTrapGradient(iChannel))
                {
                    ++nGradTrap;
                }
            }

            // scanning mdh boundaries / m_maxMdhids, part 1 (here we only update current counters)
            if (pBlock->isLabel())
            {
                ++nLBL;
                labelStateAndBookkeeping.updateLabelValues(pBlock);
                // dump_labels(arbitrarySBB.GetmaxMdhids(),"m_maxMdhids (global max boundary), in prep()");
            }
            // Count ADC events
            if (pBlock->isADC())
            {
                ++nADC;
                // if (m_nAdcLength == 0)
                //{
                //    m_nAdcLength       = pBlock->GetADCEvent().numSamples;
                //    m_nAdcDwellTime_ns = pBlock->GetADCEvent().dwellTime;
                //}
                // scanning mdh boundaries / m_maxMdhids, part 2 (here we actually update max)
                if (nLBL)
                {
                    labelStateAndBookkeeping.updateBookkeepingRecordsADC();
                }
            }

            if (pBlock->isRotation())
            {
                ++nRot;
            }
            if (pBlock->isSoftDelay())
            {
                ++nSoDe;
                /* int nNumID = pBlock->GetSoftDelayEvent().numID;
                if (nNumID >= 0 && nNumID < MAX_SOFT_DELAYS)
                {
                    if (m_vsSoDeHint.size() <= nNumID)
                    {
                        m_vsSoDeHint.resize(nNumID + 1);
                        for (int i = m_vdSoDeMin.size(); i <= nNumID; ++i)
                        {
                            m_vdSoDeMin.push_back(-1e5); // maybe there are better ways, but in the majority of cases the
                                                         // difference is not relevant
                            m_vdSoDeMax.push_back(1e5);
                            m_vdSoDeDef.push_back(0.0);
                        }
                    }
                    m_vsSoDeHint[nNumID] = pBlock->GetSoftDelayEvent().hint;
                    double dTest = double(-pBlock->GetSoftDelayEvent().offset) * pBlock->GetSoftDelayEvent().factor *
                1e-3; if (pBlock->GetSoftDelayEvent().factor >= 0)
                    {
                        if (m_vdSoDeMin[nNumID] < dTest)
                            m_vdSoDeMin[nNumID] = dTest;
                    }
                    else
                    {
                        if (m_vdSoDeMax[nNumID] > dTest)
                            m_vdSoDeMax[nNumID] = dTest;
                    }
                    // calculate default soft delay value
                    dTest = double(pBlock->GetStoredDuration() - pBlock->GetSoftDelayEvent().offset)
                            * pBlock->GetSoftDelayEvent().factor * 1e-3;
                    if (m_vdSoDeDef[nNumID] == 0.0 && dTest > 0.0)
                        m_vdSoDeDef[nNumID]
                            = dTest; // we only update the default once based on the first positive encounter
                    //
                    ExternalSequence::print_msg(
                        NORMAL_MSG,
                        std::ostringstream().flush()
                            << "block " << lI << " has soft delay " << pBlock->GetSoftDelayEvent().hint << " with ID "
                            << pBlock->GetSoftDelayEvent().numID);
                    ExternalSequence::print_msg(
                        DEBUG_MEDIUM_LEVEL,
                        std::ostringstream().flush() << "min:" << m_vdSoDeMin[nNumID] << " max:" << m_vdSoDeMax[nNumID]
                                                     << " def:" << m_vdSoDeDef[nNumID]);

                    // account for the measure time increase due to the soft delay
                    switch (arbitrarySBB.getCurrentLabelValue(ONCE))
                    {
                        case 0:
                            if (m_vmapOffsetTimeFactor_Inner.size() <= nNumID)
                                m_vmapOffsetTimeFactor_Inner.resize(nNumID + 1);
                            m_vmapOffsetTimeFactor_Inner[nNumID]
                                                        [pBlock->GetSoftDelayEvent().offset
                                                         * pBlock->GetSoftDelayEvent().factor]
                                += 1.0 / pBlock->GetSoftDelayEvent().factor;
                            break;
                        case 1:
                            if (m_vmapOffsetTimeFactor_WarmUp.size() <= nNumID)
                                m_vmapOffsetTimeFactor_WarmUp.resize(nNumID + 1);
                            m_vmapOffsetTimeFactor_WarmUp[nNumID]
                                                         [pBlock->GetSoftDelayEvent().offset
                                                          * pBlock->GetSoftDelayEvent().factor]
                                += 1.0 / pBlock->GetSoftDelayEvent().factor;
                            break;
                        case 2:
                            if (m_vmapOffsetTimeFactor_CoolDown.size() <= nNumID)
                                m_vmapOffsetTimeFactor_CoolDown.resize(nNumID + 1);
                            m_vmapOffsetTimeFactor_CoolDown[nNumID]
                                                           [pBlock->GetSoftDelayEvent().offset
                                                            * pBlock->GetSoftDelayEvent().factor]
                                += 1.0 / pBlock->GetSoftDelayEvent().factor;
                            break;
                    }
                }*/
            }

            // ExternalSequence::print_msg(DEBUG_HIGH_LEVEL, std::ostringstream().flush() <<
            //	"}*** m_llMeasureTimeInner    = " << m_llMeasureTimeInner <<
            //	"\n                 ***m_llMeasureTimeWarmUp   = " << m_llMeasureTimeWarmUp <<
            //	"\n                 ***m_llMeasureTimeCoolDown = " << m_llMeasureTimeCoolDown);

            // m_dMeasureTimeUsec += (double)pBlock->GetDuration();
            /*if (!pBlock->isSoftDelay()) // contributions of the soft delays to the measurement time are accounted for
                                        // differently, namely using m_vmapOffsetTimeFactor_Inner,
            m_vmapOffsetTimeFactor_...
            {
                switch (arbitrarySBB.getCurrentLabelValue(ONCE))
                {
                    case 0:
                        m_llMeasureTimeInner += pBlock->GetDuration();
                        break;
                    case 1:
                        m_llMeasureTimeWarmUp += pBlock->GetDuration();
                        break;
                    case 2:
                        m_llMeasureTimeCoolDown += pBlock->GetDuration();
                        break;
                }
            }*/
            llTotalDuration += int(0.5+pBlock->GetDuration());

            // free up the memory
            delete pBlock;
        }

        std::cout << "The sequence contains " << SequenceData.GetNumberOfBlocks() << " blocks with the total of\n";
        if (nRF)
            std::cout << "\t " << nRF << " RF pulses\n";
        if (nGradTrap)
            std::cout << "\t " << nGradTrap << " trapezoid gradients\n";
        if (nGradExtTrap)
            std::cout << "\t " << nGradExtTrap << " extended trapezoid gradients\n";
        if (nGradArb)
            std::cout << "\t " << nGradArb << " shaped gradients\n";
        if (nADC)
            std::cout << "\t " << nADC << " ADC events\n";
        if (nLBL)
            std::cout << "\t " << nLBL << " label extension events\n";
        if (nRot)
            std::cout << "\t " << nRot << " rotation extension events\n";
        if (nSoDe)
            std::cout << "\t " << nSoDe << " soft delay extension events\n";

        std::cout << std::endl;
        std::cout << "True total sequence duration: " << llTotalDuration << " us\n\n";
        // TODOs:
        // print info about RF pulses and arbitrary gradients, etc...

        if (labelStateAndBookkeeping.anyLabelsInUse())
        {
            std::cout << "This sequence uses LABEL extension, dumping the final results...\n";
            labelStateAndBookkeeping.dump();
        }
    }
    //std::cout << "The sequence objest should have been unloaded from memory at this point, if you see this message, it did work out.\n";
    //std::this_thread::sleep_for(5s);

	return 0;
}