//-------------------------------------------------------------------------------------------
void ReconstructShowers_TC()
{
    Log(2, "ShowRec.cpp", "--- void ReconstructShowers_TC() ---");

    //-----------------------------------------------------------------
    // Main function for reconstruction of TC "Track(s) (in) Cone" Algorithm
    // Connects tracks to ConteTube Reconstructed Object like in the CA or OI Alg.
    //-----------------------------------------------------------------

    //-----------------------------------
    // For each InitiatorBT this is
    // divided in several small parts:
    //
    // 1) Make local_gAli with cut parameters, Make GetPID of InBT and corresponding of plates
    // 2a) Loop over (whole) local_gAli, check BT for Cuts (ConeTube)
    // 2b) Loop over (whole) tracks from a linked_tracks file, try to attach these tracks to BT of these within the cone
    // 2c) Loop over (whole) attached basetracks, start ConeTube Reconstruction for these BT itsself.
    // 4) Calculate pur/eff/NBT numbers
    // 5) Fill Trees
    //-----------------------------------

    // Define Helper Variables:
    EdbPVRec* local_gAli;
    EdbSegP* InBT;
    EdbSegP* seg;
    Float_t local_gAli_pat_interim_halfsize=0;

    GLOBAL_InBTArrayEntries=GLOBAL_InBTArray->GetEntries();
    GLOBAL_INBTSHOWERNR=0;

    Float_t gAliZMax=TMath::Max(GLOBAL_gAli->GetPattern(GLOBAL_gAli->Npatterns()-1)->Z(),GLOBAL_gAli->GetPattern(0)->Z());

    if (gEDBDEBUGLEVEL>2) {
        cout << "--- --- --- Printing GLOBAL_InBTArray  MaximumZ Position " << gAliZMax << "  --- --- ---"<<endl;
        PrintShowerObjectArray(GLOBAL_InBTArray);
    }

    // Define TObjArray which contains Bastracks Clones from the tracking.
    EdbSegP * segment=0;
    EdbSegP * s2=0;
    EdbTrackP  *t  = 0;
    TFile * fil = new TFile("linked_tracks.root");
    TTree* tr= (TTree*)fil->Get("tracks");
    TClonesArray *segClonesArray= new TClonesArray("EdbSegP",60);
    int nentr = int(tr->GetEntries());
    int nseg,n0,npl;
    cout << "Of linked_tracks we have " << nentr << "  entries Total"<<endl;
    tr->SetBranchAddress("t.", &t  );
    tr->SetBranchAddress("s", &segClonesArray  );
    tr->SetBranchAddress("nseg", &nseg  );
    tr->SetBranchAddress("n0", &n0  );
    tr->SetBranchAddress("npl", &npl  );
    if (gEDBDEBUGLEVEL>2) tr->Show(0);

    cout << "Try Creating Arrays"<< endl;
    TObjArray* LOCAL_TrackBTArray= new TObjArray(nentr);
    TObjArray* LOCAL_NewAddedBTArray= new TObjArray(nentr);
    cout << "Creating Arrays done."<< endl;

    for (int i=0; i<nentr; i++ ) {
        //if (i%10!=0) continue;
        tr->GetEntry(i);

        if (gEDBDEBUGLEVEL>3) cout << i << "    t.eX Y Z " << t->X() << "   "<< t->Y() << "   "<< t->Z() << endl;

        // // Take only Basetracks from tracks which pass 3 or more plates:
        // if (npl<3) continue;


        // Take only tracks which are inside the gAli Volume:
        // Cause Linking is done with all 57 plates, but reconstruction can be done of coure on less:
        s2=(EdbSegP*) segClonesArray->At(nseg-1);
        if (s2->Z()>gAliZMax) {
            cout << " TrackEND is outside of gALi Volume!  Dont take track." << endl;
            continue;
        }

        //       CUTPARAMETER[6]=CUT_TRACKATTACH_NTRACKSEG
        //       CUT_TRACKATTACH_NTRACKSEG ranges from [0..nseg-2]
        //       nseg ranges from [2..nseg-1]

        int takeN=TMath::Min(nseg-1,int(CUT_PARAMETER[6])+1);
        //cout << " takeN =  " << takeN << endl;

        for (int hh=0; hh<takeN; hh++) {
            // Take first BT of the tracks:
            s2=(EdbSegP*) segClonesArray->At(hh);
            if (gEDBDEBUGLEVEL>3) s2->PrintNice();

            // Note: here we take all tracks, regardingless if they are in our desired range. (FP,MP,LP)
            // This we check further down:
            // Clone segment, otherwise it will be overwritten by the next.
            segment=(EdbSegP*)s2->Clone();
            // Add it to LOCAL_TrackBTArray
            LOCAL_TrackBTArray->Add(segment);
        }
        //cout << "for (int hh=0; hh<takeN; hh++) done."<<endl;

    }
    delete tr;
    delete fil;

    if (gEDBDEBUGLEVEL>1) cout << "--- Filled " << LOCAL_TrackBTArray->GetEntries()  << " Segments into LOCAL_TrackBTArray."<<endl;


    //-----------------------------------------------------------------
    // Since GLOBAL_InBTArray is filled in ascending ordering by zpositon
    // We use the descending loop to begin with BT with lowest z first.
    // This part of the Algorithm is the same as OI Alg part:
    for (Int_t i=GLOBAL_InBTArrayEntries-1; i>=0; --i) {

        //-----------------------------------
        // CounterOutPut
        if (gEDBDEBUGLEVEL==2) if ((i%1)==0) cout << GLOBAL_InBTArrayEntries <<" InBT in total, still to do:"<<Form("%4d",i)<< "\r\r\r\r"<<flush;
        //-----------------------------------

        //-----------------------------------
        // Get InitiatorBT from GLOBAL_InBTArray
        InBT=(EdbSegP*)GLOBAL_InBTArray->At(i);
        //--------
        GLOBAL_InBT_E=InBT->P();
        GLOBAL_InBT_TanTheta=TMath::Sqrt(InBT->TX()*InBT->TX()+InBT->TY()*InBT->TY());
        GLOBAL_InBT_Flag=InBT->Flag();
        GLOBAL_InBT_MC=InBT->MCEvt();
        //--------
        Int_t local_NBT=0;
        Int_t local_NBTMC=0;
        Int_t local_NBTallMC=0;
        Int_t local_NBTeMC=0;
        float_t local_pure=-1;
        float_t local_purall=-1;
        Int_t npat_int=0;
        Int_t npat_total=0;
        Int_t npatN=0;
        Int_t npat_Neff=0;
        Int_t NBT_Neff=0;
        Int_t NBTMC_Neff=0;
        Int_t NBTMCe_Neff=0;
        //--------

        if (gEDBDEBUGLEVEL>2) {
            cout << endl << endl << "--- Starting Shower for Number " << i << " now: "<<endl;
            InBT->PrintNice();
            cout << " InBT->MCEvt()  "<< InBT->MCEvt()  << "   InBT->P()  " << InBT->P()  << endl;
        }
        //-----------------------------------

        //-----------------------------------
        // 1) Make local_gAli with cut parameters:
        //-----------------------------------
        local_gAli = TransformEdbPVRec(GLOBAL_gAli, InBT);
        //local_gAli = TransformEdbPVRec_SA(GLOBAL_gAli, InBT); // Especially for SA-Algorithm to set the size of the local_gAli also.

        // Add InBT to GLOBAL_ShowerSegArray
        GLOBAL_ShowerSegArray -> Add(InBT);
        //-----------------------------------


        //-----------------------------------
        // 2a) Loop over (whole) local_gAli, check BT for Cuts
        //-----------------------------------
        Int_t local_gAli_npat=local_gAli->Npatterns();
        if (gEDBDEBUGLEVEL>2) cout << "--- local_gAli_npat=  " << local_gAli_npat << endl;

        // Loop over all plates of local_gAli, since this is already
        // extracted with the right numbers of plates...
        for (Int_t patterloop_cnt=local_gAli_npat-1; patterloop_cnt>=0; --patterloop_cnt) {
            if (gEDBDEBUGLEVEL>3) cout << "--- --- Doing patterloop_cnt= " << patterloop_cnt << endl;

            for (Int_t btloop_cnt=0; btloop_cnt<local_gAli->GetPattern(patterloop_cnt)->GetN(); ++btloop_cnt) {
                seg = (EdbSegP*)local_gAli->GetPattern(patterloop_cnt)->GetSegment(btloop_cnt);
                if (gEDBDEBUGLEVEL>3) seg->PrintNice();

                // Now apply cut conditions: TrackCone Alg  --------------------
                // (for first iteration, thi is equivalent to OI Alg Implementation:
                //         if (GetdeltaRWithPropagation(seg,InBT)>100) continue;
                //         if (TMath::Abs(seg->Z()-InBT->Z())>5000) continue;
                if (!GetConeOrTubeDistanceToInBT(seg, InBT, CUT_PARAMETER[0], CUT_PARAMETER[1])) continue;
                if (!FindPrecedingBTsSingleThetaAngle(seg, InBT, local_gAli, GLOBAL_ShowerSegArray)) continue;
                // end of    cut conditions: TrackCone IMPLEMENTATION Alg  --------------------

                // If we arrive here, Basetrack seg has passed criteria
                // and is then added to the shower array:
                // Check if its not the InBT which is already added:
                AddBTToArrayWithCeck(seg,  GLOBAL_ShowerSegArray);
            }

        }  // end of loop over all plates of local_gAli
//     if (gEDBDEBUGLEVEL>2) PrintShowerObjectArray(GLOBAL_ShowerSegArray);


        Float_t NBT_before=0;
        Float_t NBT_after=0;
        Float_t NBTMC_before=0;
        Float_t NBTMC_after=0;
        Float_t pur_before=0;
        Float_t pur_after=0;

        // Quickly calculate NEW NBT and pur------------------------
        EdbSegP* interim=0;
        for (int count=0; count<GLOBAL_ShowerSegArray->GetEntries(); count++) {
            NBT_before++;
            interim=(EdbSegP*)GLOBAL_ShowerSegArray->At(count);
            if (interim->MCEvt()>0) NBTMC_before++;
        }
        pur_before=(Float_t)NBTMC_before/(Float_t)NBT_before;
        // Quickly calculate NEW NBT and pur------------------------



        //-----------------------------------
        // 2b) Loop over (whole) tracks from a linked_tracks file, try to attach these tracks to BT of these within the cone
        // wARNING,wARNING: the linked tracks root file BaseTrack PIDs  and the PIDs of the BTs in the eAli DO NOT
        // NECESSARILY have to be identical, so please, make comparisons rather on Z instead on PID.
        //-----------------------------------


        //cout << "=============================================================================="<<endl;
        //cout << "=" << endl;
        //cout << "=   NOW try to attach tracks to the reconstructed showers:" << endl;

        // Clear LOCAL_NewAddedBTArray to have only same MCEvts in it:
        LOCAL_NewAddedBTArray->Clear();

        int nentr_LOCAL_TrackBTArray=LOCAL_TrackBTArray->GetEntries();
        int nentr_GLOBAL_ShowerSegArray=GLOBAL_ShowerSegArray->GetEntries();

        int  nentr2=nentr_LOCAL_TrackBTArray;

        EdbSegP* tryAttachedSegment=0;
        EdbSegP* alreadyTakenSegment=0;

        for (Int_t icount=0; icount<nentr2; icount++ ) {
            //for(int icount=nentr2/2; icount<(nentr2/2)+2; icount++ ) {

            tryAttachedSegment=(EdbSegP*)LOCAL_TrackBTArray->At(icount);

            // But take only segemnst with same MC Number as InBT of course!!
            if (tryAttachedSegment->MCEvt()>0 && tryAttachedSegment->MCEvt()!=InBT->MCEvt()) continue;

            //cout << "Try to attac now segment (if array) " << icount << " with: " << endl;
            //tryAttachedSegment->PrintNice();
            //cout << "to any of the already attached segments: " << endl;

            for (Int_t j=0; j<GLOBAL_ShowerSegArray->GetEntries(); j++ ) {
                alreadyTakenSegment=(EdbSegP*)GLOBAL_ShowerSegArray->At(j);

                Float_t dRTEST= GetdeltaRWithPropagation(tryAttachedSegment,alreadyTakenSegment);
                Float_t dMinDist=GetMinimumDist(tryAttachedSegment,alreadyTakenSegment);
                Float_t dMinDist_DT=GetdeltaThetaSingleAngles(tryAttachedSegment,alreadyTakenSegment);

                /*
                if (gEDBDEBUGLEVEL>2) {
                cout << "     Test GLOBAL_ShowerSegArray segment (if array) " << j << " with: " << endl;
                alreadyTakenSegment->PrintNice();
                cout << "     to connect: deltaRWith Propagation (tryAttachedSegment,alreadyTakenSegment) : ";
                cout << dRTEST;;
                cout << "     GetMinimumDist(tryAttachedSegment,alreadyTakenSegment) : " << dMinDist << " " << dMinDist_DT << endl;
                }
                */

                // The testing trackBT has to satisfy both criteria:
                // minimum Distance and deltaThetaSingleAngels
                // CUT_PARAMETER[4]=CUT_TRACKATTACH_DISTMIN
                // CUT_PARAMETER[5]=CUT_TRACKATTACH_DTAN_MAX
                if (dMinDist< CUT_PARAMETER[4] && dMinDist_DT< CUT_PARAMETER[5]) {
                    LOCAL_NewAddedBTArray->Add(tryAttachedSegment);
                    //if (gEDBDEBUGLEVEL>2) cout << "Found a close connection. Add tryAttachedSegment to  LOCAL_NewAddedBTArray"<<endl;
                    break;
                }

            } // of (int j=0;)
        } // of (int icount=0; );

        //cout << "LOCAL_NewAddedBTArray->GetEntries()  = " << LOCAL_NewAddedBTArray->GetEntries() << endl;
        //cout << "=" << endl;
        //cout << "=============================================================================="<<endl;


        //-----------------------------------
        // 2c) Loop over (whole) local_gAli, again, but now with ONLY the new LOCAL_NewAddedBTArray as starting BTs
        // 2c) and only if the local_gAli Z position id greater equal than LOCAL_NewAddedBTArray z position:
        // 2c) Do this for every new added BaseTrack:
        // 2c) Rember that not the BaseTracks of the linked_tracks.root file are NOT given into ShowerOject array,
        // 2c) Only the intrinsic ones which come out from the gAli directly.
        // 2c) This is to avoid duplication of tracks, since the may not be totally equal in X,Y,Z
        // 2c) For example if tracking was done with small different par sets...
        //-----------------------------------
        //cout << " // 2c) Loop over (whole) local_gAli, again, but now with ONLY the new LOCAL_NewAddedBTArray as starting BTs"<<endl;
        //cout << " // 2c) LOCAL_NewAddedBTArray->GetEntries()  " << LOCAL_NewAddedBTArray->GetEntries() << endl;
        //PrintShowerObjectArray(LOCAL_NewAddedBTArray);

        for (Int_t icount=0; icount<LOCAL_NewAddedBTArray->GetEntries(); icount++) {

            tryAttachedSegment=(EdbSegP*)LOCAL_NewAddedBTArray->At(icount);
            if (tryAttachedSegment->MCEvt()>0 && tryAttachedSegment->MCEvt()!=InBT->MCEvt()) continue;

            //cout << "// 2c) Loop over (whole) local_gAli, again, but now with ONLY the new LOCAL_NewAddedBTArray as starting BTs"<< endl;
            //cout << "// 2c) and only if the local_gAli Z position id greater equal than LOCAL_NewAddedBTArray z position:"<< endl;
            //cout << " Doing this for  tryAttachedSegment  " << tryAttachedSegment->X() << " " << tryAttachedSegment->Y() << " " << tryAttachedSegment->Z() << " " << tryAttachedSegment->TX() << " " <<  tryAttachedSegment->TY() << " " << tryAttachedSegment->MCEvt() << " " << endl;

            local_gAli_npat=local_gAli->Npatterns();
            // cout << "--- local_gAli_npat=  " << local_gAli_npat << endl;

            // Loop over all plates of local_gAli, since this is already
            // extracted with the right numbers of plates...
            for (Int_t patterloop_cnt=local_gAli_npat-1; patterloop_cnt>=0; --patterloop_cnt) {

                float_t local_gAliZ=local_gAli->GetPattern(patterloop_cnt)->Z();

                //         cout << "--- --- Doing patterloop_cnt= " << patterloop_cnt << endl;
                //         cout << "--- --- local_gAli->GetPattern(patterloop_cnt)->Z()= " << local_gAliZ << endl;
                //         cout << "--- --- Is tryAttachedSegment->Z() gtreater equal local_gAliZ()= " << tryAttachedSegment->Z() << " " << local_gAliZ << endl;

                // We do not want to search for local_gAli which has lower Z than the attacet BT itself
                //if (tryAttachedSegment->Z() < local_gAliZ) continue;
                //if (tryAttachedSegment->Z() > local_gAliZ) continue;

                int matches_found=0;
                //cout << "--- --- YES its greater, searching BTs of localgAli now for matches  ?? ? GREATER OR SMALLER ???"<<endl;


                for (Int_t btloop_cnt=0; btloop_cnt<local_gAli->GetPattern(patterloop_cnt)->GetN(); ++btloop_cnt) {
                    seg = (EdbSegP*)local_gAli->GetPattern(patterloop_cnt)->GetSegment(btloop_cnt);
                    if (gEDBDEBUGLEVEL>3) seg->PrintNice();

                    // Now what was formerly InBT is now the new tryAttachedSegment:
                    // cout << " comparing now Segment "  << seg << "  (Z= " << seg->Z() << ")  with the tryAttachedSegment (Z= " << tryAttachedSegment->Z() << ") GetdeltaRWithPropagation: "<< GetdeltaRWithPropagation(seg,tryAttachedSegment) << endl;

                    // Now apply cut conditions: TrackCone Alg  --------------------
                    if (!GetConeOrTubeDistanceToInBT(seg, tryAttachedSegment, CUT_PARAMETER[0], CUT_PARAMETER[1])) continue;
                    // cout << "Passed GetConeOrTubeDistanceToInBT"<< endl;
                    if (!FindPrecedingBTsSingleThetaAngleTCDEBUG(seg, tryAttachedSegment, local_gAli, GLOBAL_ShowerSegArray)) continue;
                    // cout << "Passed FindPrecedingBTsSingleThetaAngleTCDEBUG"<< endl;
                    // end of    cut conditions: TrackCone IMPLEMENTATION Alg  --------------------

                    // If we arrive here, Basetrack seg has passed criteria
                    // and is then added to the shower array:
                    // Check if its not the InBT which is already added:
                    Bool_t isBTadded=kFALSE;
                    isBTadded = AddBTToArrayWithCeck(seg,GLOBAL_ShowerSegArray);

                    if (gEDBDEBUGLEVEL>2) {
                        if (isBTadded) cout << "...Yes this was a new one, and is added to shower array! ..."<< endl;
                        if (isBTadded) seg->PrintNice();
                    }
                    if (isBTadded) matches_found++;

                } // of  (Int_t btloop_cnt=0; );

                //if (matches_found>0) cout << "--- --- In this plate we have found _ADDITIONAL_ (more than zero) BTs: "<< matches_found << endl;

                // Calculate Efficency numbers only for the last element of the loop,
                // otherwise we add to much numbers...
                if (icount==LOCAL_NewAddedBTArray->GetEntries()-1) {
                    // Calc BT density around shower:
                    EdbPattern* pat_interim=local_gAli->GetPattern(patterloop_cnt);
                    CalcTrackDensity(pat_interim,local_gAli_pat_interim_halfsize,npat_int,npat_total,npatN);
                    // Calc TrackNumbers for plate for efficency numbers:
                    CalcEfficencyNumbers(pat_interim, InBT->MCEvt(), NBT_Neff, NBTMC_Neff,NBTMCe_Neff);
                }

            }  // for (Int_t patterloop_cnt=local_gAli_npat-1; patterloop_cnt>=0; --patterloop_cnt)  // end of loop over all plates of local_gAli

        } // for(int icount=0; icount<OCAL_NewAddedBTArray->GetEntries(); icount++ )
        LOCAL_NewAddedBTArray->Clear();

//      if (gEDBDEBUGLEVEL>2) PrintShowerObjectArray(GLOBAL_ShowerSegArray);


        // Quickly calculate NEW NBT and pur------------------------
        NBT_after=0;
        NBTMC_after=0;
        for (int count=0; count<GLOBAL_ShowerSegArray->GetEntries(); count++) {
            NBT_after++;
            interim=(EdbSegP*)GLOBAL_ShowerSegArray->At(count);
            if (interim->MCEvt()>0) NBTMC_after++;
        }
        pur_after=(Float_t)NBTMC_after/(Float_t)NBT_after;
        // Quickly calculate NEW NBT and pur------------------------

        if (gEDBDEBUGLEVEL>2) {
            cout << "  Difference between before and after: " << endl;
            cout << "  NBT =  " << NBT_after <<  "   Before: "<< NBT_before << endl;
            cout << "  NBTMC =  " << NBTMC_after <<  "   Before: "<< NBTMC_before << endl;
            cout << "  pur =  " << pur_after <<  "   Before: "<< pur_before << endl;
        }

        //-----------------------------------
        // 4) Calculate pur/eff/NBT numbers,
        // not needed when only reconstruction
        // done:
        //-----------------------------------
        if (cmd_OUTPUTLEVEL>=2 || cmd_OUTPUTLEVEL==0 ) {
            Int_t NBT=0;
            Int_t NBTMC=0;
            Int_t NBTallMC=0;
            Int_t NBTeMC=0;
            Double_t  eff, purall, pure;
            CalcEffPurOfShower2(GLOBAL_ShowerSegArray, NBT, NBTMC, NBTallMC, NBTeMC, purall, pure, NBT_Neff, NBTMC_Neff,NBTMCe_Neff);

            // Fill only for MC Event:
            if (GLOBAL_InBT_MC>0) {
                GLOBAL_EvtBT_Flag=GLOBAL_EvtBT_FlagArray[GLOBAL_InBT_MC];
                GLOBAL_EvtBT_MC=GLOBAL_EvtBT_MCArray[GLOBAL_InBT_MC];
                GLOBAL_EvtBT_E=GLOBAL_EvtBT_EArray[GLOBAL_InBT_MC];
                GLOBAL_EvtBT_TanTheta=GLOBAL_EvtBT_TanThetaArray[GLOBAL_InBT_MC];
            }
            GLOBAL_trckdens=shower_trackdensb;
        }


        //-----------------------------------
        // 5) Fill Tree:
        //-----------------------------------
        TREE_ShowRecEff->Fill();
        if (gEDBDEBUGLEVEL>3) TREE_ShowRecEff->Show(TREE_ShowRecEff->GetEntries()-1);


        //-----------------------------------
        // 6a) Transfer ShowerArray to treebranchTreeEntry:
        //-----------------------------------
        if (cmd_OUTPUTLEVEL>0) {
            TransferShowerObjectArrayIntoEntryOfTreebranchShowerTree(TREE_ShowShower,GLOBAL_ShowerSegArray);
        }


        //------------------------------------
        // Reset and delete important things:
        // also  to avoid memory problems ...
        //-----------------------------------

        GLOBAL_ShowerSegArray->Clear();
        if (gEDBDEBUGLEVEL>3) cout << "--- ---GLOBAL_ShowerSegArray->GetEntries(): "<< GLOBAL_ShowerSegArray->GetEntries() << endl;
        delete local_gAli;
        local_gAli=0;
        ++GLOBAL_INBTSHOWERNR;
        //------------------------------------
    }
    // end of loop over GLOBAL_InBTArrayEntries
    //-----------------------------------------------------------------


    delete LOCAL_TrackBTArray;
    delete LOCAL_NewAddedBTArray;

    if (gEDBDEBUGLEVEL==2) cout << endl<<flush;
    if (gEDBDEBUGLEVEL>3) cout << "---TREE_ShowRecEff->GetEntries() ... " << TREE_ShowRecEff->GetEntries() << endl;
    if (gEDBDEBUGLEVEL>3) cout << "---GLOBAL_INBTSHOWERNR ... " << GLOBAL_INBTSHOWERNR<< endl;


    return;
}
//-------------------------------------------------------------------------------------------
