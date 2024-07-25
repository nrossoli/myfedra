//-------------------------------------------------------------------------------------------
void ReconstructShowers_CL()
{
    Log(2, "ShowRec.cpp", "--- void ReconstructShowers_CL() ---");

    //-----------------------------------------------------------------
    // Main function for reconstruction of "CLuster" Algorithm
    //-----------------------------------------------------------------
    //
    // For each InitiatorBT this is
    // divided in several small parts:
    // 1) Make local_gAli with cut parameters, Make GetPID of InBT and corresponding of plates
    // 3) Loop over (whole) local_gAli, check BT for Cuts
    // 4) Calculate pur/eff/NBT numbers
    // 5) Fill Trees
    //-----------------------------------

    //-----------------------------------
    //   0)  Predefinitions
    //-----------------------------------

    // Set ClusterBoxSize of X,Y in microns: (f.e. 50x50)
    Float_t BOXSIZE=50;
    Float_t ClusterBoxSizeX=BOXSIZE;
    Float_t ClusterBoxSizeY=ClusterBoxSizeX;
    Float_t ClusterBoxSize[]= {ClusterBoxSizeX,ClusterBoxSizeY};

    // Set GLOBAL Boundaries of the BrickSize, depends on your MC Simulation and on your scanned BG size
    // Later we will set the frame size, which can be f.e. 1mmX1mm, 3mmX3mm ...
    // This does not need to be changed...
    Float_t LeftX=-11000;
    Float_t RightX=21000;
    Float_t LeftY=-11000;
    Float_t RightY=16000;

    cout << " ......."<<endl;
    // Set SquareSize of The gAliSub -> Coorseponds to the totral scan Area:
    Float_t VIEWSIZE=1000;
    Float_t gAliSub_Length_X=VIEWSIZE;
    Float_t gAliSub_Length_Y=VIEWSIZE; //1mm

    // Calculate how many NbinsX,Y are needed:
    Int_t NbinsX=Int_t((RightX-LeftX)/ClusterBoxSizeX);
    Int_t NbinsY=Int_t((RightY-LeftY)/ClusterBoxSizeY);


    cout << " ......."<<endl;

    Int_t NPLATES=GLOBAL_gAli->Npatterns();
    Int_t NbinsZ=NPLATES;

    // Calculate how many NbinsX,Y are needed:
    Int_t NbinsX_gAliSub=Int_t((gAliSub_Length_X)/ClusterBoxSizeX);
    Int_t NbinsY_gAliSub=Int_t((gAliSub_Length_Y)/ClusterBoxSizeY);

    // Arrays for the Spectrum source and dests...
    // Float_t ** source = new float *[NbinsX];
    Double_t ** source = new double *[NbinsX];
    for (Int_t cnti=0; cnti<NbinsX; cnti++) {
        // source[cnti]=new float[NbinsY];
        source[cnti]=new double[NbinsY];
    }
    // Float_t ** dest = new float *[NbinsX];
    Double_t ** dest = new double *[NbinsX];
    for (Int_t cnti=0; cnti<NbinsX; cnti++) {
        // dest[cnti]=new float[NbinsY];
        dest[cnti]=new double[NbinsY];
    }

    // Arrays for the Spectrum source and dests...
    // Unused ??? Not to be found in any other part of the code ...
    // Comment out!
    /*
    Float_t ** source_gAliSub = new float *[NbinsX_gAliSub];
    for (Int_t cnti=0; cnti<NbinsX_gAliSub; cnti++) {
        source_gAliSub[cnti]=new float[NbinsY_gAliSub];
    }
    Float_t ** dest_gAliSub = new float *[NbinsX_gAliSub];
    for (Int_t cnti=0; cnti<NbinsX_gAliSub; cnti++) {
        dest_gAliSub[cnti]=new float[NbinsY_gAliSub];
    }
    */

    // Fill Extrapolated BT with direction into next plate?
    Bool_t FillTANData=kTRUE;
    cout << " ......."<<endl;

    // Plates to be grouped Together:
    Int_t  PLATESPERGROUP=1;
    Int_t  FIRSTPLATE=0;
    Int_t LASTPLATE=0;
    Int_t  MAXPLATE=GLOBAL_gAli->Npatterns()-1;
    Int_t NGroupedPLATES=ceil( (Double_t)NPLATES/(Double_t)PLATESPERGROUP);

    // Spectrum Position Peaks
    Float_t* fPositionX;
    Float_t* fPositionY;
    Int_t fNPeaks;

    // Start Position...
    Float_t StartPosX,StartPosY,StartPosTX,StartPosTY;

    // Calculate the z-positions for the grouped plates:
    Float_t ZPosGroupedPlates[50];
    ///=== DEBUG===
    cout << "gAli->GetPattern(firstplate)->GetSegment(0)->Z;  "<<GLOBAL_gAli->GetPattern(0)->Z()<<endl;
    cout << "gAli->GetPattern(lastplate)->GetSegment(0)->Z;  "<<GLOBAL_gAli->GetPattern(GLOBAL_gAli->Npatterns()-1)->Z()<<endl;
    ///===END DEBUG===

    for (Int_t i=0; i<NGroupedPLATES; i++) {
        ZPosGroupedPlates[i]=GLOBAL_gAli->GetPattern(0)->Z()+PLATESPERGROUP*1300*i;
        ///===DEBUG===
        cout << "i PLATESPERGROUP NGroupedPLATES ZPosGroupedPlates[i] " << i << " " <<  PLATESPERGROUP << " " << NGroupedPLATES  << " " << ZPosGroupedPlates[i] << " " << endl;
    }

    ///====== DEBUG===
    if (gEDBDEBUGLEVEL==2) {
        cout << "NbinsX = " << NbinsX << endl;
        cout << "NbinsY = " << NbinsY << endl;
        cout << "NbinsZ = " << NbinsZ << endl;
        cout << "NbinsX_gAliSub = " << NbinsX_gAliSub << endl;
        cout << "NbinsY_gAliSub = " << NbinsY_gAliSub << endl;
        cout << "FIRSTPLATE=  "<<FIRSTPLATE<<endl;
        cout << "LASTPLATE=  "<<LASTPLATE<<endl;
        cout << "MAXPLATE=  "<<MAXPLATE<<endl;
        cout << "NPLATES=  "<<NPLATES<<endl;
        cout << "NGroupedPLATES=  "<<NGroupedPLATES<<endl;
        cout << "PLATESPERGROUP=  "<<PLATESPERGROUP<<endl;
    }
    ///===END DEBUG===

    //=========================================
    // Create the 2Dim  Histograms:
    Reco_CL_BuildGlobalHistogramsOnHeap();
    // Create the 2Dim  Histograms:
    Hist2DimOnlyBGAllPlates= new TH2F();
    Hist2DimOnlySimOneEventAllPlates= new TH2F();
    Hist2DimBGAndSimOneEventAllPlates= new TH2F();
    // Create the 3Dim  Histograms:
    Hist3DimOnlyBG = new TH3F();
    Hist3DimRecoEvent_gAli = new TH3F();
    //=========================================
    // Allocate TSectrum2 pointer on the heap, and try to speed up program....
    spectrum2dim= new TSpectrum2();
    spectrum2dim->Print();
    //=========================================


    // Define Helper Variables:
    EdbPVRec* local_gAli;
    EdbSegP* InBT;
    EdbSegP* seg;
    Float_t local_gAli_pat_interim_halfsize=0;

    GLOBAL_InBTArrayEntries=GLOBAL_InBTArray->GetEntries();
    GLOBAL_INBTSHOWERNR=0;
    cout << "GLOBAL_InBTArrayEntries = " << GLOBAL_InBTArrayEntries << endl;

    TClonesArray  *as;
    EdbSegP  *sa;
    EdbPattern  *as_pattern;
    EdbPattern  *as_pattern_sub;


    //-----------------------------------------------------------------
    // Since GLOBAL_InBTArray is filled in ascending ordering by zpositon
    // We use the descending loop to begin with BT with lowest z first.
    //   for (Int_t i=GLOBAL_InBTArrayEntries-1; i>=0; --i) {
    for (Int_t i=GLOBAL_InBTArrayEntries-1; i>=GLOBAL_InBTArrayEntries-1; --i) {

        //-----------------------------------
        // CounterOutPut
        if (gEDBDEBUGLEVEL==2) if ((i%1)==0) cout << GLOBAL_InBTArrayEntries <<" InBT in total, still to do:"<<Form("%4d",i)<< "\r\r\r\r"<<flush;
        if (gEDBDEBUGLEVEL==1) {
            int modulo=GLOBAL_InBTArrayEntries/20;
            if ((i%modulo)==0) cout << i <<" : 5% more done"<<endl;
        }
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
        }
        //-----------------------------------

        //-----------------------------------
        // 1) Make local_gAli with cut parameters:
        //-----------------------------------
        local_gAli = TransformEdbPVRec(GLOBAL_gAli, InBT);
        // Add InBT to GLOBAL_ShowerSegArray
        GLOBAL_ShowerSegArray -> Add(InBT);
        //-----------------------------------

        //-----------------------------------
        // 1a) Set Some Variables
        //-----------------------------------
        StartPosX=InBT->X();
        StartPosY=InBT->Y();
        StartPosTX=InBT->TX();
        StartPosTY=InBT->TY();
        //-----------------------------------

        //-----------------------------------
        //  1b) Now SetBins for the 2D, 3D Histograms for THIS (specific) Event
        //-----------------------------------

        cout << "Reset now some other histos"<<endl;
        for (Int_t h=0; h<50; h++) {
            Hist2DimOnlyBGOneGroupedPlate[h]->Reset();
            Hist2DimOnlySimOneEventOneGroupedPlate[h]->Reset();
            Hist2DimBGAndSimOneEventOneGroupedPlate[h]->Reset();
        }

        cout << "Reset now some other SetNameTitle "<<endl;
        Hist2DimOnlyBGAllPlates->SetNameTitle("Hist2DimOnlyBGAllPlates","Hist2DimOnlyBGAllPlates");
        Hist2DimOnlySimOneEventAllPlates->SetNameTitle("Hist2DimOnlySimOneEventAllPlates","Hist2DimOnlySimOneEventAllPlates");
        Hist2DimBGAndSimOneEventAllPlates->SetNameTitle("Hist2DimBGAndSimOneEventAllPlates","Hist2DimBGAndSimOneEventAllPlates");
        Hist2DimOnlyBGAllPlates->SetBins(  NbinsX_gAliSub, StartPosX-gAliSub_Length_X/2.0, StartPosX +gAliSub_Length_X/2.0, NbinsY_gAliSub, StartPosY-gAliSub_Length_Y/2.0, StartPosY +gAliSub_Length_Y/2.0);
        Hist2DimOnlySimOneEventAllPlates->SetBins(  NbinsX_gAliSub, StartPosX-gAliSub_Length_X/2.0, StartPosX +gAliSub_Length_X/2.0, NbinsY_gAliSub, StartPosY-gAliSub_Length_Y/2.0, StartPosY +gAliSub_Length_Y/2.0);
        Hist2DimBGAndSimOneEventAllPlates->SetBins(  NbinsX_gAliSub, StartPosX-gAliSub_Length_X/2.0, StartPosX +gAliSub_Length_X/2.0, NbinsY_gAliSub, StartPosY-gAliSub_Length_Y/2.0, StartPosY +gAliSub_Length_Y/2.0);

        cout << "Reset now some other SetNameTitle "<<endl;
        Hist3DimOnlyBG->SetNameTitle("Hist3DimOnlyBG","Hist3DimOnlyBG");
        Hist3DimOnlyBG->SetBins( NbinsX_gAliSub, StartPosX -gAliSub_Length_X/2.0, StartPosX +gAliSub_Length_X/2.0, NbinsY_gAliSub, StartPosY -gAliSub_Length_Y/2.0, StartPosY +gAliSub_Length_Y/2.0, NGroupedPLATES, (FIRSTPLATE-1)*1300, NbinsZ*1300);
        Hist3DimRecoEvent_gAli->SetNameTitle("Hist3DimRecoEvent_gAli","Hist3DimRecoEvent_gAli");
        Hist3DimRecoEvent_gAli->SetBins( NbinsX, LeftX, RightX, NbinsY, LeftY, RightY, NGroupedPLATES, (FIRSTPLATE-1)*1300, (LASTPLATE)*1300);

        cout << "Reset now some other SetNameTitle "<<endl;
        Hist2DimOnlyBGAllPlates->Reset();
        Hist2DimOnlySimOneEventAllPlates->Reset();
        Hist2DimBGAndSimOneEventAllPlates->Reset();
        Hist3DimOnlyBG->Reset();
        Hist3DimRecoEvent_gAli->Reset();

        //-----------------------------------
        // 1c) Define the Histograms Bounds and Bins,
        // 1c) Since Bounds are depending on h ...
        //-----------------------------------
        Float_t InBT_X=InBT->X();
        Float_t InBT_TX=InBT->TX();
        Float_t InBT_X_Extrapolated;
        Float_t InBT_Y=InBT->Y();
        Float_t InBT_TY=InBT->TY();
        Float_t InBT_Y_Extrapolated;
        Float_t InBT_Z=InBT->Z();
        for (Int_t h=0; h<NGroupedPLATES; h++) {
            //cout << "--------------------"<<endl;
            Float_t zdiff=ZPosGroupedPlates[h]-InBT_Z;
            InBT_X_Extrapolated=InBT_X+InBT_TX*zdiff;
            InBT_Y_Extrapolated=InBT_Y+InBT_TY*zdiff;
            cout << "h InBT_X InBT_TX InBT_Z ZPosGroupedPlates[h] zdiff  InBT_X_Extrapolated " << h << " " <<  InBT_X << " " << InBT_TX << " " << InBT_Z << " " << ZPosGroupedPlates[h] << " " << zdiff << " " <<InBT_X_Extrapolated << endl;
            cout << "h InBT_Y InBT_TY InBT_Z ZPosGroupedPlates[h] zdiff  InBT_Y_Extrapolated " << h << " " <<  InBT_Y << " " << InBT_TY << " " << InBT_Z << " " << ZPosGroupedPlates[h] << " " << zdiff << " " <<InBT_Y_Extrapolated << endl;
            Hist2DimOnlyBGOneGroupedPlate[h]->SetBins(NbinsX_gAliSub, InBT_X_Extrapolated-gAliSub_Length_X/2.0,InBT_X_Extrapolated+gAliSub_Length_X/2.0, NbinsY_gAliSub,InBT_Y_Extrapolated-gAliSub_Length_Y/2.0,InBT_Y_Extrapolated+gAliSub_Length_Y/2.0);
            Hist2DimOnlySimOneEventOneGroupedPlate[h]->SetBins(NbinsX_gAliSub, InBT_X_Extrapolated-gAliSub_Length_X/2.0,InBT_X_Extrapolated+gAliSub_Length_X/2.0, NbinsY_gAliSub,InBT_Y_Extrapolated-gAliSub_Length_Y/2.0,InBT_Y_Extrapolated+gAliSub_Length_Y/2.0);
            Hist2DimBGAndSimOneEventOneGroupedPlate[h]->SetBins(NbinsX_gAliSub, InBT_X_Extrapolated-gAliSub_Length_X/2.0,InBT_X_Extrapolated+gAliSub_Length_X/2.0, NbinsY_gAliSub,InBT_Y_Extrapolated-gAliSub_Length_Y/2.0,InBT_Y_Extrapolated+gAliSub_Length_Y/2.0);
        }
        //=========================================



        //-----------------------------------
        // 1c) Loop over (whole) local_gAli, Fill the Histograms  (from firstplate up to lastplate)
        //-----------------------------------
        Int_t val, val_tan;
        Int_t local_gAli_npat=local_gAli->Npatterns();
        if (gEDBDEBUGLEVEL>2) cout << "--- local_gAli_npat=  " << local_gAli_npat << endl;

        // Loop over all plates of local_gAli, since this is already
        // extracted with the right numbers of plates...
        for (Int_t patterloop_cnt=local_gAli_npat-1; patterloop_cnt>=0; --patterloop_cnt) {
            if (gEDBDEBUGLEVEL>3) cout << "--- --- Doing patterloop_cnt= " << patterloop_cnt << endl;

            as = (TClonesArray*)local_gAli->GetPattern(patterloop_cnt)->GetSegments();
            as_pattern= (EdbPattern*)local_gAli->GetPattern(patterloop_cnt);

            cout << "===DEBUG   patterloop_cnt=  " << patterloop_cnt << "   Z pos (k) " << local_gAli->GetPattern(patterloop_cnt)->Z() << endl;  // get Z pos of the plate k
            cout << "===DEBUG   (patterloop_cnt) EdbSegments in the Pattern (total): " << as->GetEntries() << endl;

            for (Int_t h=0; h<as->GetEntries(); h++) {
                sa = (EdbSegP*)( as->At(h) );
                ///if (sa->X()<LeftX || sa->X()>RightX || sa->Y()<LeftY || sa->Y()>RightY) continue;  /// DEBUG
                val=Reco_CL_AssignZValueToGroup(sa->Z(), ZPosGroupedPlates[0], NGroupedPLATES, PLATESPERGROUP);
                cout << "entry, val  "<<h << "  "  << val<<endl;
                Hist2DimBGAndSimOneEventOneGroupedPlate[val]->Fill(sa->X(),sa->Y());
                if (sa->MCEvt()<0) Hist2DimOnlyBGOneGroupedPlate[val]->Fill(sa->X(),sa->Y());
                if (FillTANData) {
                    val_tan=Reco_CL_AssignZValueToGroup(sa->Z()+1300, ZPosGroupedPlates[0], NGroupedPLATES, PLATESPERGROUP);
                    cout << "(TAN) entry, val_tan  "<<h << "  "  << val_tan<<endl;
                    if (val_tan>=NGroupedPLATES) {
                        /// --val_tan; // Project the last on beack to its owns
                        ++val_tan; // Project the last on beack to its owns
                        cout << "(TAN) val_tan changed to "<<h << "  "  << val_tan<<endl;
                    }
                    Hist2DimBGAndSimOneEventOneGroupedPlate[val_tan]->Fill(sa->X()+sa->TX()*1300,sa->Y()+sa->TY()*1300);
                    if (sa->MCEvt()<0) Hist2DimOnlyBGOneGroupedPlate[val_tan]->Fill(sa->X()+sa->TX()*1300,sa->Y()+sa->TY()*1300);
                }
                if (sa->MCEvt()<0) continue;
                if (sa->MCEvt()!=GLOBAL_InBT_MC) continue;
                cout << "===   ==DEBUG   Filling i with MCEvt  "<< i <<"("<< GLOBAL_InBT_MC <<")  X,Y,Z: "<<sa->X()<<" "<<sa->Y()<<" "<<sa->Z()<<"  "<<sa->MCEvt()<< "   to val= " <<val <<endl;
                Hist2DimOnlySimOneEventOneGroupedPlate[val]->Fill(sa->X(),sa->Y());
                if (FillTANData) {
                    Hist2DimOnlySimOneEventOneGroupedPlate[val_tan]->Fill(sa->X()+sa->TX()*1300,sa->Y()+sa->TY()*1300);
                }
            }

        } //of Loop over all plates of local_gAli,


        cout<<endl<<endl;
        cout <<"Loop over the Grouped plates and search for spectrum peaks in each groupedPlate."<<endl;
        ///-------------------------------------------
        /// ====   Now loop over the Grouped plates and search for spectrum peaks in each groupedPlate...
        for (Int_t h=0; h<NGroupedPLATES; h++) {

            /// Set the bin arrays for the Spectrum...
            for (Int_t cnti = 0; cnti < NbinsX_gAliSub; cnti++) {
                for (Int_t cntj = 0; cntj < NbinsY_gAliSub; cntj++) {
                    source[cnti][cntj] = 0; // Reset source before filling it new...
                    dest[cnti][cntj] = 0; // Reset dest before filling it new...
                    source[cnti][cntj] = Hist2DimBGAndSimOneEventOneGroupedPlate[h]->GetBinContent(cnti + 1,cntj + 1);
                }
            }
            ///--------------------------------------------
            /// Fit each grouped Plate with Spectrum...
            Int_t nfound=0;
            cout << " Do now peak search..."<<endl;
            nfound = spectrum2dim->SearchHighRes(source, dest, NbinsX_gAliSub, NbinsY_gAliSub, 2, 20, kTRUE, 10, kFALSE, 5);// TO BE OPTIMIZED...
            cout << " Peak search finished. Go on."<<endl;
            ///--------------------------------------------
            spectrum_interim=(TH2F*)Hist2DimBGAndSimOneEventOneGroupedPlate[h]->Clone();
            spectrum_interim->Reset();
            ///.---------------
            Int_t interimsbin=0;
            TAxis *xax;
            TAxis *yax;
            TAxis *zax;
            Double_t xax_bin_value;
            Double_t yax_bin_value;
            Double_t zax_bin_value;
            Double_t value;

            /// Fill interimSpectrum For Drawing... without the THRESHOLD_SMOOTHED_DEST Cut to find Maximum
            for (Int_t ii = 0; ii < NbinsX_gAliSub; ii++) {
                for (Int_t jj = 0; jj < NbinsY_gAliSub; jj++) {
                    value=dest[ii][jj];
                    spectrum_interim->SetBinContent(ii + 1,jj + 1, value);
                }
            }
            Float_t THRESHOLD_SMOOTHED_DEST=0.5;
            cout << "THRESHOLD_SMOOTHED_DEST  spectrum_interim->GetMaximum()  : " << THRESHOLD_SMOOTHED_DEST << "  " << spectrum_interim->GetMaximum()<< endl;
            cout << "----------_"<<endl;

            Double_t spec_int_maximum=spectrum_interim->GetMaximum();

            /// Fill interimSpectrum For Drawing...only when Entry is > THRESHOLD_SMOOTHED_DEST * Maximum
            /// THRESHOLD_SMOOTHED_DEST is in % Units of Maxium...
            for (Int_t ii = 0; ii < NbinsX_gAliSub; ii++) {
                for (Int_t jj = 0; jj < NbinsY_gAliSub; jj++) {
                    // cout << " h ii jj THRESHOLD_SMOOTHED_DEST dest[ii][jj] "<< h << " " <<  ii << "  " << jj << "  " << THRESHOLD_SMOOTHED_DEST  << "  " << dest[ii][jj] << endl;
                    if (dest[ii][jj]<THRESHOLD_SMOOTHED_DEST*spec_int_maximum) continue;
                    value=dest[ii][jj];
                    spectrum_interim->SetBinContent(ii + 1,jj + 1, value);
                    xax_bin_value=spectrum_interim->GetXaxis()->GetBinCenter(ii + 1);
                    yax_bin_value=spectrum_interim->GetYaxis()->GetBinCenter(jj + 1);
                    zax_bin_value=ZPosGroupedPlates[h];
                    //      cout <<" xax_bin   ii jj (xyz) value  "<<  xax_bin_value  << "  " <<  ii <<  "  " << jj<< "  " <<xax_bin_value  << "  " << yax_bin_value<< "  " << zax_bin_value<<"    " << value << endl;
                    Hist3DimRecoEvent_gAli->Fill(xax_bin_value,yax_bin_value,zax_bin_value, value);
                }
            }
            cout << "THRESHOLD_SMOOTHED_DEST  spectrum_interim->GetMaximum()  : " << THRESHOLD_SMOOTHED_DEST << "  " << spectrum_interim->GetMaximum()<< endl;
            cout << "----------_"<<endl;

            //AllGroupsSimPlusBGittedSpectrum->cd(h+1);
            //spectrum_interim->DrawCopy("colz");
        }
        /// end of loop over the Grouped plates
        ///-------------------------------------------





        return;

        //-----------------------------------
        // 2) Loop over (whole) local_gAli, check BT for Cuts
        //-----------------------------------
        //     Int_t local_gAli_npat=local_gAli->Npatterns();
        if (gEDBDEBUGLEVEL>2) cout << "--- local_gAli_npat=  " << local_gAli_npat << endl;

        // Loop over all plates of local_gAli, since this is already
        // extracted with the right numbers of plates...
        for (Int_t patterloop_cnt=local_gAli_npat-1; patterloop_cnt>=0; --patterloop_cnt) {
            if (gEDBDEBUGLEVEL>3) cout << "--- --- Doing patterloop_cnt= " << patterloop_cnt << endl;

            for (Int_t btloop_cnt=0; btloop_cnt<local_gAli->GetPattern(patterloop_cnt)->GetN(); ++btloop_cnt) {
                seg = (EdbSegP*)local_gAli->GetPattern(patterloop_cnt)->GetSegment(btloop_cnt);
                if (gEDBDEBUGLEVEL>3) seg->PrintNice();

                // Now apply cut conditions: CT ConeTube Alg  --------------------
                if (!GetConeOrTubeDistanceToInBT(seg, InBT, CUT_PARAMETER[0], CUT_PARAMETER[1])) continue;
                if (!FindPrecedingBTs(seg, InBT, local_gAli, GLOBAL_ShowerSegArray)) continue;
                // end of    cut conditions: CT ConeTube Alg  --------------------

                // If we arrive here, Basetrack seg has passed criteria
                // and is then added to the shower array:
                // Check if its not the InBT which is already added:
                if (seg->X()==InBT->X()&&seg->Y()==InBT->Y()) {
                    ;    // do nothing;
                }
                else {
                    GLOBAL_ShowerSegArray -> Add(seg);
                }
            }
            // Calc BT density around shower:
            EdbPattern* pat_interim=local_gAli->GetPattern(patterloop_cnt);
            CalcTrackDensity(pat_interim,local_gAli_pat_interim_halfsize,npat_int,npat_total,npatN);

            // Calc TrackNumbers for plate for efficency numbers:
            CalcEfficencyNumbers(pat_interim, InBT->MCEvt(), NBT_Neff, NBTMC_Neff,NBTMCe_Neff);
        }
        // end of loop over all plates of local_gAli
        if (gEDBDEBUGLEVEL>2) PrintShowerObjectArray(GLOBAL_ShowerSegArray);

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





Int_t Reco_CL_AssignZValueToGroup(Double_t z, Double_t z0, Int_t NGroupedPLATES, Int_t PLATESPerGroup)
{
    // z:  Position of BT to be checked,
    // z0: ZStart if for example FIRSTPLATE z=3900 (PLATE 4) instead of 0
    Int_t actplate=(Int_t)(z-z0)/1300+1;
    Int_t actgroup=(actplate-1)/PLATESPerGroup;
    //  cout << "z  z0 actplate "<< z << "   " << z0 << "  " <<actplate << endl;
    //  cout << "actplate  actgroup "<< actplate << "   " << actgroup << endl;
    //This Assignment seems ok.
    return actgroup;
}

void Reco_CL_BuildGlobalHistogramsOnHeap()
{
    TString HistoNameTitle="XXX";
    for (Int_t h=0; h<50; h++) {
        //    cout << "Creating    Hist2DimOnlyBG_Groupe_ " << h << endl;
        HistoNameTitle=TString(Form("Hist2DimOnlyBG_Groupe_%d",h));
        Hist2DimOnlyBGOneGroupedPlate[h]=new TH2F(HistoNameTitle,HistoNameTitle,10,0,1,10,0,1);
        //    cout << "Creating    Hist2DimOnlySimOneEvent_Groupe_ " << h << endl;
        HistoNameTitle=TString(Form("Hist2DimOnlySimOneEvent_Groupe_%d",h));
        Hist2DimOnlySimOneEventOneGroupedPlate[h]=new  TH2F(HistoNameTitle,HistoNameTitle,10,0,1,10,0,1);
        //    cout << "Creating    Hist2DimBGAndSimOneEvent_Groupe_ " << h << endl;
        HistoNameTitle=TString(Form("Hist2DimBGAndSimOneEvent_Groupe_%d",h));
        Hist2DimBGAndSimOneEventOneGroupedPlate[h]=new  TH2F(HistoNameTitle,HistoNameTitle,10,0,1,10,0,1);
    }
    return;
}

