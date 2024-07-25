//-------------------------------------------------------------------------------------------
void ReconstructShowers_BW()
{
    Log(2, "ShowRec.cpp", "--- void ReconstructShowers_BW() ---");

    //-----------------------------------------------------------------
    // Main function for reconstruction of "BackWard" Algorithm
    //-----------------------------------------------------------------

    //-----------------------------------
    // For each InitiatorBT this is
    // divided in several small parts:
    //
    // 1) Make local_gAli with cut parameters, Make GetPID of InBT and corresponding of plates
    // 3) Loop over (whole) local_gAli, check BT for Cuts
    // 4) Calculate pur/eff/NBT numbers
    // 5) Fill Trees
    //-----------------------------------

    //-----------------------------------
    // Algorithm Iteration Steps:
    //
    //  0) Start Adjust Shower Axis, Virtual Vertex
    //  0) Calc dR/dT/dIP(trk,trk) to same plate.
    //  ---Loop up to firstplate;
    //  //  ---Loop for deltaN (plate)=1,2,3:
    // //  //  1) Calc dR/dT/dIP(trk,trk) to backward plate.
    // //  //  1) Check for already in shower array; add BT to shower
    // //  //  1) Adjust Shower Axis, Adjust Virtual Vertex
    //  //  ---Loop for deltaN (plate)=1,2,3:
    //  ---Loop up to firstplate;
    //  2) Adjust Shower Axis, Virtual Vertex.
    //
    //-----------------------------------


    // Define Helper Variables:
    EdbPVRec* local_gAli;
    EdbSegP* InBT;
    EdbSegP* seg;
    EdbSegP* segShower;
    EdbSegP* ShowerAxis;
    EdbVertex* VirtualVertex;
    Float_t local_gAli_pat_interim_halfsize=0;
    Double_t dR,dT,dminDist;

    GLOBAL_InBTArrayEntries=GLOBAL_InBTArray->GetEntries();
    GLOBAL_INBTSHOWERNR=0;

    //-----------------------------------------------------------------
    // Since GLOBAL_InBTArray is filled in ascending ordering by zpositon
    // We use the descending loop to begin with BT with lowest z first.
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
        Int_t show_NP=0;
        Int_t npat_Neff=0;
        Int_t NBT_Neff=0;
        Int_t NBTMC_Neff=0;
        Int_t NBTMCe_Neff=0;
        //--------

        if (gEDBDEBUGLEVEL!=2) {
            cout << endl << endl << "--- Starting Shower for Number " << i << " (MCEvt=" << InBT->MCEvt() << ") now: "<<endl;
            InBT->PrintNice();
            cout << "InBT->P() " << InBT->P() << endl;
        }
        //-----------------------------------

        //-----------------------------------
        // 1) Make local_gAli with cut parameters:
        //-----------------------------------
        local_gAli = TransformEdbPVRec_BackWard(GLOBAL_gAli, InBT);
        // Add InBT to GLOBAL_ShowerSegArray
        GLOBAL_ShowerSegArray -> Add(InBT);
        if (gEDBDEBUGLEVEL>2) cout << "--- TransformEdbPVRec done."<<endl;
        //-----------------------------------


        //-----------------------------------
        // 2) Start at plate for Initator BT:
        //-----------------------------------
        Int_t local_gAli_npat=local_gAli->Npatterns();
        Int_t btloop_cnt_N=0;
        //   if (gEDBDEBUGLEVEL>2)
        cout << "--- Starting BACKWARD Reco... ---" << endl;
        if (gEDBDEBUGLEVEL>2) cout << "--- local_gAli_npat=  " << local_gAli_npat << endl;

        Int_t patterloop_cnt_InBT=-1;
        Float_t patterloop_cnt_InBT_Z=0;
        for (Int_t patterloop_cnt=0; patterloop_cnt<local_gAli_npat; patterloop_cnt++) {
            if (gEDBDEBUGLEVEL>3) cout << "--- --- Doing patterloop_cnt= " << patterloop_cnt << endl;
            // if pattern is after the Initiator BT->Z() position we stop:
            if (local_gAli->GetPattern(patterloop_cnt)->Z()<InBT->Z()) continue;
            patterloop_cnt_InBT=patterloop_cnt;
            patterloop_cnt_InBT_Z=local_gAli->GetPattern(patterloop_cnt)->Z();
        }
        if (gEDBDEBUGLEVEL>2) cout << "--- patterloop_cnt_InBT=  " << patterloop_cnt_InBT << endl;
        if (gEDBDEBUGLEVEL>2) cout << "--- patterloop_cnt_InBT_Z=  " << patterloop_cnt_InBT_Z << endl;

        // Check wich Steptype is needed for BackwardStepping:
        // Convention: -1 if Patter(0)_Z<Pattern(npat)_Z
        // Convention: +1 if Patter(0)_Z>Pattern(npat)_Z
        Int_t StepNr=0;
        if (gEDBDEBUGLEVEL>2) cout << "local_gAli->GetPattern(0)->Z()  " <<  local_gAli->GetPattern(0)->Z() << endl;
        if (gEDBDEBUGLEVEL>2) cout << "local_gAli->GetPattern(local_gAli_npat-1)->Z()  " <<  local_gAli->GetPattern(local_gAli_npat-1)->Z() << endl;
        if (local_gAli->GetPattern(0)->Z() > local_gAli->GetPattern(local_gAli_npat-1)->Z()) {
            StepNr=1;
        }
        else {
            StepNr=-1;
        }
        if (gEDBDEBUGLEVEL>2) cout << "--- StepNr for going backwards (decreasing Z) to next plate=  " << StepNr << endl;

        // return;
        Bool_t FP_reached=kFALSE;
        Int_t patterloop_cnt=patterloop_cnt_InBT;
        if (StepNr==1)  patterloop_cnt=patterloop_cnt_InBT;
        if (StepNr==-1)  cout << "DONT KNOW HERE WHICH    patterloop_cnt=local_gAli_npat-1;   to calculate!!! " << endl;

        // Loop over all plates of local_gAli, since this is already
        // extracted with the right numbers of plates...
        while (!FP_reached) {
            if (gEDBDEBUGLEVEL>2) cout << " FP NOT YET REACHED .... << doing patterloop_cnt = " << patterloop_cnt << endl;
            btloop_cnt_N=local_gAli->GetPattern(patterloop_cnt)->GetN();
            if (gEDBDEBUGLEVEL>2) cout << " FP NOT YET REACHED .... << with  btloop_cnt_N = " << btloop_cnt_N << endl;

            // Loop over all BTs in the actual plate. Calculate dR/dT/dIP to
            // all BTs which are already in the shower (at beginning it is
            // only the InBT.
            for (Int_t btloop_cnt=0; btloop_cnt<btloop_cnt_N; ++btloop_cnt) {
                seg = (EdbSegP*)local_gAli->GetPattern(patterloop_cnt)->GetSegment(btloop_cnt);
                if (gEDBDEBUGLEVEL>3) seg->PrintNice();

                Bool_t add_seg=kFALSE;
                if (gEDBDEBUGLEVEL>3) cout << "btloop_cnt (ot of ) " << btloop_cnt << "  ( "  << btloop_cnt_N << " ) " << endl;

                // Now apply cut conditions: BACKWARD Alg  --------------------
                // Check for dTheta to tracks already in shower:
                Int_t actualEntries=GLOBAL_ShowerSegArray->GetEntries();
                if (gEDBDEBUGLEVEL>3) cout << "actualEntries of  GLOBAL_ShowerSegArray:  " << actualEntries << endl;
                for (int i=0; i<actualEntries; ++i) {
                    segShower=(EdbSegP*)GLOBAL_ShowerSegArray->At(i);
                    Float_t dZ=TMath::Abs(segShower->Z()-seg->Z());
                    dT=GetdeltaThetaSingleAngles(seg,segShower); // DO NOT TAKE GetdeltaTheta since this is calculation based on absolute theta differences (not on TX,TY relative ones) !!
                    dR=GetdR(seg,segShower);
                    dminDist=GetMinimumDist(seg,segShower);
                    if (gEDBDEBUGLEVEL>3) cout << "btloop_cnt i  dT  dR  dminDist dZ: " << btloop_cnt << "  " << i << "  " << dT << "  " << dR << "  " << dminDist << "  " << dZ << endl;
                    if (TMath::Abs(dZ)>CUT_PARAMETER[3]) continue;
                    if (TMath::Abs(dT)>CUT_PARAMETER[1]) continue;
                    if (TMath::Abs(dR)>CUT_PARAMETER[2]) continue;
                    if (TMath::Abs(dminDist)>CUT_PARAMETER[0]) continue;
                    if (gEDBDEBUGLEVEL>3) cout << "try to add this BT (if not already in there...) "<< endl;
                    add_seg=kTRUE;
                    break;
                } // of for (int i=0; i<actualEntries; ++i)
                if (add_seg) {
                    AddBTToArrayWithCeck(seg,GLOBAL_ShowerSegArray);
                    EdbSegP* BT_1=(EdbSegP*)GLOBAL_ShowerSegArray->At(0);
                    EdbSegP* BT_2=(EdbSegP*)GLOBAL_ShowerSegArray->At(GLOBAL_ShowerSegArray->GetEntries()-1);
                    show_NP=TMath::Max(show_NP,TMath::Abs(BT_1->PID()-BT_2->PID())+1);
                    //cout << "show_NP = " << show_NP << endl;
                }
            } // of for (Int_t btloop_cnt=0; btloop_cnt<btloop_cnt_N; ++btloop_cnt)

            // Now goto next plate:
            patterloop_cnt=patterloop_cnt+StepNr;
            if (patterloop_cnt<0) FP_reached=kTRUE;
            if (patterloop_cnt>=local_gAli_npat) FP_reached=kTRUE;
            // Or we stop also if the number of plates is more than cmd_NP
            if (show_NP>=cmd_NP) FP_reached=kTRUE;
        } // while (!FP_reached)


        if (gEDBDEBUGLEVEL>2) PrintShowerObjectArray(GLOBAL_ShowerSegArray);
        if (cmd_OUTPUTLEVEL>=2 || cmd_OUTPUTLEVEL==0 ) {
            Int_t NBT=0;
            Int_t NBTMC=0;
            Int_t NBTallMC=0;
            Int_t NBTeMC=0;
            Double_t  purall, pure;
            CalcEffPurOfShower(GLOBAL_ShowerSegArray, NBT, NBTMC, NBTallMC, NBTeMC, purall, pure);
            GLOBAL_EvtBT_Flag=GLOBAL_EvtBT_FlagArray[i];
            GLOBAL_EvtBT_MC=GLOBAL_EvtBT_MCArray[i];
            GLOBAL_EvtBT_E=GLOBAL_EvtBT_EArray[i];
            GLOBAL_EvtBT_TanTheta=GLOBAL_EvtBT_TanThetaArray[i];
            GLOBAL_EvtBT_Flag=GLOBAL_EvtBT_FlagArray[i];
        }

        if (gEDBDEBUGLEVEL>2) cout << " Check if its sorted: GLOBAL_ShowerSegArray->At(0)->Z() should be the lowest one...continue;" << endl;
        if (gEDBDEBUGLEVEL>2) cout << " IsShowerSortedZ(GLOBAL_ShowerSegArray)  " << IsShowerSortedZ(GLOBAL_ShowerSegArray) << endl;
        SortShowerZ(GLOBAL_ShowerSegArray);

        // After Sort Shower we can build axis right now...
        ShowerAxis=BuildShowerAxis(GLOBAL_ShowerSegArray);
        ShowerAxis->PrintNice();
        Float_t mindTTT=999999;
        Int_t mindTTT_h=0;
        for (int h=0; h<GLOBAL_ShowerSegArray->GetEntries(); ++h) {
            segShower=(EdbSegP*)GLOBAL_ShowerSegArray->At(h);
            //   segShower->PrintNice();
            if (segShower->Z()!=ShowerAxis->Z()) continue;
            Float_t dTTT = GetdeltaTheta(ShowerAxis,segShower); //cout << "dTTT  = " << dTTT << endl;
            if (dTTT<mindTTT) {
                mindTTT=dTTT;
                mindTTT_h=h;
            }
        }

        // Now here comes the forward (==downstream==inbeamdirection) reconstruction
        if (gEDBDEBUGLEVEL>2) cout << "--- Starting FORWARD Reco... ---" << endl;

        // Just take the BT which is closest to shower axis:
        ShowerAxis=(EdbSegP*)GLOBAL_ShowerSegArray->At(mindTTT_h);
        GLOBAL_ShowerSegArray->Clear();
        // The First we have to add by hand because in the reco routine there has to be at least one BT to check.
        GLOBAL_ShowerSegArray->Add(ShowerAxis);
        if (gEDBDEBUGLEVEL>3) ShowerAxis->PrintNice();

        for (Int_t patterloop_cnt=0; patterloop_cnt<local_gAli_npat; patterloop_cnt++) {
            if (gEDBDEBUGLEVEL>3) cout << "--- --- Doing patterloop_cnt= " << patterloop_cnt << endl;
            // if pattern is after the ShowerAxis BT->Z() position we stop:
            if (local_gAli->GetPattern(patterloop_cnt)->Z()<ShowerAxis->Z()) continue;
            patterloop_cnt_InBT=patterloop_cnt;
            patterloop_cnt_InBT_Z=local_gAli->GetPattern(patterloop_cnt)->Z();
        }
        if (gEDBDEBUGLEVEL>2) cout << "--- patterloop_cnt_InBT=  " << patterloop_cnt_InBT << endl;
        if (gEDBDEBUGLEVEL>2) cout << "--- patterloop_cnt_InBT_Z=  " << patterloop_cnt_InBT_Z << endl;

        // InvertSpeNr now for forward Step:
        Int_t StepNrForward=StepNr*-1;
        StepNr=StepNrForward;
        show_NP=0;

        Bool_t LP_reached=kFALSE;
        patterloop_cnt=patterloop_cnt_InBT;

        while (!LP_reached) {
            if (gEDBDEBUGLEVEL>2) cout << " LP NOT YET REACHED .... << doing patterloop_cnt = " << patterloop_cnt << endl;
            btloop_cnt_N=local_gAli->GetPattern(patterloop_cnt)->GetN();
            if (gEDBDEBUGLEVEL>2) cout << " LP NOT YET REACHED .... << with  btloop_cnt_N = " << btloop_cnt_N << endl;

            // Loop over all BTs in the actual plate. Calculate dR/dT/dIP to
            // all BTs which are already in the shower (at beginning it is
            // only the InBT.
            for (Int_t btloop_cnt=0; btloop_cnt<btloop_cnt_N; ++btloop_cnt) {
                seg = (EdbSegP*)local_gAli->GetPattern(patterloop_cnt)->GetSegment(btloop_cnt);
                if (gEDBDEBUGLEVEL>3) seg->PrintNice();
                Bool_t add_seg=kFALSE;
                if (gEDBDEBUGLEVEL>3) cout << "btloop_cnt (ot of ) " << btloop_cnt << "  ( "  << btloop_cnt_N << " ) " << endl;

                // Now apply cut conditions: BACKWARD Alg  --------------------
                // Check for dTheta to tracks already in shower:
                Int_t actualEntries=GLOBAL_ShowerSegArray->GetEntries();
                if (gEDBDEBUGLEVEL>3) cout << "actualEntries of  GLOBAL_ShowerSegArray:  " << actualEntries << endl;
                for (int i=0; i<actualEntries; ++i) {
                    segShower=(EdbSegP*)GLOBAL_ShowerSegArray->At(i);
                    Float_t dZ=TMath::Abs(segShower->Z()-seg->Z());
                    dT=GetdeltaThetaSingleAngles(seg,segShower); // DO NOT TAKE GetdeltaTheta since this is calculation based on absolute theta differences (not on TX,TY relative ones) !!
                    dR=GetdR(seg,segShower);
                    dminDist=GetMinimumDist(seg,segShower);
                    if (gEDBDEBUGLEVEL>3) cout << "btloop_cnt i  dT  dR  dminDist dZ: " << btloop_cnt << "  " << i << "  " << dT << "  " << dR << "  " << dminDist << "  " << dZ << endl;
                    if (TMath::Abs(dZ)>CUT_PARAMETER[7]) continue;
                    if (TMath::Abs(dT)>CUT_PARAMETER[5]) continue;
                    if (TMath::Abs(dR)>CUT_PARAMETER[6]) continue;
                    if (TMath::Abs(dminDist)>CUT_PARAMETER[4]) continue;
                    if (gEDBDEBUGLEVEL>3) cout << "try to add this BT (if not already in there...) "<< endl;
                    add_seg=kTRUE;
                    break;
                } // of for (int i=0; i<actualEntries; ++i)
                if (add_seg) {
                    AddBTToArrayWithCeck(seg,GLOBAL_ShowerSegArray);
                    EdbSegP* BT_1=(EdbSegP*)GLOBAL_ShowerSegArray->At(0);
                    EdbSegP* BT_2=(EdbSegP*)GLOBAL_ShowerSegArray->At(GLOBAL_ShowerSegArray->GetEntries()-1);
                    show_NP=TMath::Max(show_NP,TMath::Abs(BT_1->PID()-BT_2->PID())+1);
                    //cout << "show_NP = " << show_NP << endl;
                }
            } // of for (Int_t btloop_cnt=0; btloop_cnt<btloop_cnt_N; ++btloop_cnt)


            // Calc BT density around shower:
            EdbPattern* pat_interim=local_gAli->GetPattern(patterloop_cnt);
            CalcTrackDensity(pat_interim,local_gAli_pat_interim_halfsize,npat_int,npat_total,npatN);

            // Calc TrackNumbers for plate for efficency numbers:
            CalcEfficencyNumbers(pat_interim, InBT->MCEvt(), NBT_Neff, NBTMC_Neff,NBTMCe_Neff);

            // Now goto next plate:
            patterloop_cnt=patterloop_cnt+StepNr;
            if (patterloop_cnt<0) LP_reached=kTRUE;
            if (patterloop_cnt>=local_gAli_npat) LP_reached=kTRUE;
            // Or we stop also if the number of plates is more than cmd_NP
            if (show_NP>=cmd_NP) LP_reached=kTRUE;

        } // of while (!LP_reached)
        ///   return;  --------------------------------------------------------------------------

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

    if (gEDBDEBUGLEVEL==2) cout << endl<<flush;
    if (gEDBDEBUGLEVEL>3) cout << "---TREE_ShowRecEff->GetEntries() ... " << TREE_ShowRecEff->GetEntries() << endl;
    if (gEDBDEBUGLEVEL>3) cout << "---GLOBAL_INBTSHOWERNR ... " << GLOBAL_INBTSHOWERNR<< endl;


    return;
}
//-------------------------------------------------------------------------------------------
