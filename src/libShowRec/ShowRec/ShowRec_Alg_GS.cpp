//-------------------------------------------------------------------------------------------
void ReconstructShowers_GS()
{
    Log(2, "ShowRec.cpp", "--- void ReconstructShowers_GS() ---");
    //-----------------------------------------------------------------
    // Main function for reconstruction of "Gamma Search" Algorithm
    //-----------------------------------------------------------------

    //-----------------------------------
    // For each InitiatorBT this is
    // divided in several small parts:
    //
    // 1) Make local_gAli with cut parameters, Make GetPID of InBT and corresponding of plates
    // 3) Loop over (whole) local_gAli, search and check BT pairs for suts
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

    Bool_t IsFirstLoopCount=kTRUE;
    Int_t LastGlobalMCEventNr=1;

    EdbVertex* vtx=new EdbVertex();

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
        // 0) Set vtx according to MC event,
        //    or to a half plate backward propagated segment.
        //    Do this only if the value -GBMC is NOT set,
        //    (in the case -GBMC is set, we should know the vertex anyway).
        //-----------------------------------
        //
        if (GLOBAL_InBT_MC>LastGlobalMCEventNr) LastGlobalMCEventNr=GLOBAL_InBT_MC;
        //
        //
        if (GLOBAL_InBT_MC<0&& cmd_GBMC<0) {
            // This is an "equalization" effect. If we dont put this in, we will have BG base
            // tracks always extrapolated to themselves, but SIM basetracks extrapolated to
            // the vertex. This leads to a wrong deltaZ distribution which spoils the
            // ANN Training.
            if( IsFirstLoopCount==kTRUE) cout << "Option: GLOBAL_InBT_MC<0&& cmd_GBMC<0"<<endl;
            vtx->SetXYZ(GLOBAL_VtxArrayX[LastGlobalMCEventNr],GLOBAL_VtxArrayY[LastGlobalMCEventNr],GLOBAL_VtxArrayZ[LastGlobalMCEventNr]);
        }
        //
        else if (GLOBAL_InBT_MC<0 && cmd_GBMC==0) {
            if( IsFirstLoopCount==kTRUE) cout << "Option: GLOBAL_InBT_MC<0 && cmd_GBMC==0"<<endl;
            vtx->SetXYZ(InBT->X()-650*InBT->TX(),InBT->Y()-650*InBT->TY(),InBT->Z()-650);
        }
        else if (GLOBAL_InBT_MC<0 && cmd_GBMC>0) {
            if( IsFirstLoopCount==kTRUE)  cout << "Option: GLOBAL_InBT_MC<0 && cmd_GBMC>0"<<endl;
            vtx->SetXYZ(GLOBAL_VtxArrayX[cmd_GBMC],GLOBAL_VtxArrayY[cmd_GBMC],GLOBAL_VtxArrayZ[cmd_GBMC]);
        }
        else if (GLOBAL_InBT_MC>0 && GLOBAL_IsBrickTreePGunInfo==kFALSE) {
            if( IsFirstLoopCount==kTRUE) cout << "Option: GLOBAL_InBT_MC>0 && GLOBAL_IsBrickTreePGunInfo==kFALSE"<<endl;
            vtx->SetXYZ(InBT->X()-650*InBT->TX(),InBT->Y()-650*InBT->TY(),InBT->Z()-650);
        }
        else {
            if( IsFirstLoopCount==kTRUE) cout << "Option: else"<<endl;
            vtx->SetXYZ(GLOBAL_VtxArrayX[GLOBAL_InBT_MC],GLOBAL_VtxArrayY[GLOBAL_InBT_MC],GLOBAL_VtxArrayZ[GLOBAL_InBT_MC]);
        }
        vtx->SetMC(GLOBAL_InBT_MC);

        if (gEDBDEBUGLEVEL>2) {
            cout << "The vtx info for this MC event: X:Y:Z:MC:  " << vtx->X()  << " " << vtx->Y() << " " << vtx->Z() << " " << GLOBAL_InBT_MC << endl;
            cout << "The vtx info for this MC event: IP(INBT,vtx):  " << CalcIP(InBT,vtx)<< endl;
        }
        //-----------------------------------



        //-----------------------------------
        // 1) Make local_gAli with cut parameters:
        //-----------------------------------
        local_gAli = TransformEdbPVRec(GLOBAL_gAli, InBT);

        // IN THIS ALGORITHM WE DO NOT YET AUTOMATICALLY ADD THE FIRST BT
        // SINCE WE NEED AN ADDIITONAL VERTEX CUT OF THIS FIRST BT TO CHECK!
        // // Add InBT to GLOBAL_ShowerSegArray
        // GLOBAL_ShowerSegArray -> Add(InBT);
        //-----------------------------------


        //-----------------------------------
        // 2) Loop over (whole) local_gAli, check InitiatorBT
        // 2) compatible with a segment forming a e+e- pair:
        // 2) (Loop over all plates of local_gAli, since this is already
        // 2) extracted with the right numbers of plates...)
        //-----------------------------------
        Int_t local_gAli_npat=local_gAli->Npatterns();
        Int_t btloop_cnt_N=0;

        ///============================================================================================
        ///========================  CODE FROM EdbShowAlg_GS  FROM libShowRec =========================
        ///============================================================================================



        Int_t npat=local_gAli->Npatterns();
        Int_t pat_one_bt_cnt_max,pat_two_bt_cnt_max=0;
        EdbPattern* pat_one=0;
        EdbPattern* pat_two=0;
        EdbSegP* Segment=0;
        EdbSegP* Segment2=0;
        Float_t distZ,IP_Pair_To_InBT,IP_Pair_To_InBT_SegSum;
        Float_t IP_Pair_To_InBT_Seg2;


        Segment = InBT;
        IP_Pair_To_InBT=CalcIP(Segment, vtx);
        /// Change with respect to libShowRec: here we assume that Segment will always be
        /// the Initiator BaseTrack and Segment2 is the other segment to check.







        // Loop over pattern for the first BT of the pairs:
        // Start first with the pattern with the lowest Z position.
        pat_one=local_gAli->GetPatternZLowestHighest(1);
        Float_t pat_one_Z=pat_one->Z();
        pat_one_bt_cnt_max=pat_one->GetN();

        for (Int_t pat_one_cnt=0; pat_one_cnt<npat; ++pat_one_cnt) {

            if (pat_one_cnt>0) {
                pat_one=(EdbPattern*)local_gAli->NextPattern(pat_one_Z,1);
                pat_one_Z=pat_one->Z();
                pat_one_bt_cnt_max=pat_one->GetN();
            }

            // Check if pattern Z() equals the InBZ->Z(), since we wanna have the
            // pattern one the pattern to contain the InBT:
            distZ=pat_one->Z()-InBT->Z();
            if (TMath::Abs(distZ)>10) continue;
            //cout << "distZ (pat_one->Z()-InBT->Z())= " << distZ << endl;

            // Check if InBT fulfills criteria for IP to vertex:
            if (IP_Pair_To_InBT>CUT_PARAMETER[0]) continue;
            // Check if InBT fulfills criteria for z-diff to vertex:
            if ((InBT->Z()-vtx->Z())>CUT_PARAMETER[3]) continue;

            // Now here we can add InBT since it passed also the vertex cut.
            // Therefore, the Reconstructed Shower has always InBT as first BT stored.
            if (GLOBAL_ShowerSegArray->GetEntries()==0) {
                GLOBAL_ShowerSegArray->Add(InBT);
                if (gEDBDEBUGLEVEL>2) {
                    cout << "I have added the first InBT  " << InBT << "  to GLOBAL_ShowerSegArray." << endl;
                    InBT->PrintNice();
                }
            }

            // Check if pattern dist Z to Vtx  is ok:
            distZ=pat_one->Z()-vtx->Z();
            // Z distance has to be greater zero, cause the InBT
            // and other pair BTs shall come downstream the vertex:
            if (distZ<0) continue;
            //cout << "distZ (pat_one->Z()-vtx->Z();) = " << distZ << endl;

            if (gEDBDEBUGLEVEL>2) cout << "Searching patterns: pat_one_cnt=" << pat_one_cnt << "  pat_one->Z() = " << pat_one->Z() << " pat_one_bt_cnt_max= "<< pat_one_bt_cnt_max <<endl;


            // Loop over pattern for the second BT of the pairs:
            //
            //cout << "// Loop over pattern for the second BT of the pairs: "<< endl;

            pat_two=local_gAli->GetPatternZLowestHighest(1);
            Float_t pat_two_Z=pat_two->Z();
            pat_two_bt_cnt_max=pat_two->GetN();

            for (Int_t pat_two_cnt=0; pat_two_cnt<npat; ++pat_two_cnt) {

                if (pat_two_cnt>0) {
                    pat_two=(EdbPattern*)local_gAli->NextPattern(pat_two_Z,1);
                    pat_two_Z=pat_two->Z();
                    pat_two_bt_cnt_max=pat_two->GetN();
                }

                // PID diff of two plates may be maximum [0..PidDIFFN]
                if (TMath::Abs(pat_one_cnt-pat_two_cnt)>CUT_PARAMETER[5]) continue;

                // pattern two should come downstream pattern one:
                if (pat_two->Z()<pat_one->Z()) continue;


                if (gEDBDEBUGLEVEL>2) cout << " Searching patterns: pat_two_cnt=" << pat_two_cnt << "  pat_two->Z() = " << pat_two->Z() << " pat_two_bt_cnt_max= "<< pat_two_bt_cnt_max <<endl;


                for (Int_t pat_one_bt_cnt=0; pat_one_bt_cnt<pat_one_bt_cnt_max; ++pat_one_bt_cnt) {
                    /// Segment =  (EdbSegP*)pat_one->GetSegment(pat_one_bt_cnt);
                    /// Segment = InBT;
                    /// Change with respect to libShowRec: here we assume that Segment will always be
                    /// the Initiator BaseTrack.

                    for (Int_t pat_two_bt_cnt=0; pat_two_bt_cnt<pat_two_bt_cnt_max; ++pat_two_bt_cnt) {
                        Segment2 = (EdbSegP*)pat_two->GetSegment(pat_two_bt_cnt);

                        // Ceck if segments are not (by chance) the same:
                        if (Segment2==Segment) continue;
                        if (Segment2->ID()==Segment->ID()&&Segment2->PID()==Segment->PID()) continue;
                        if (IsSameSegment(Segment2,Segment)) continue;



                        /// At first:  Check for already duplicated pairings:
                        /// if (CheckPairDuplications(Segment->PID(),Segment->ID(),Segment2->PID(),Segment2->ID(), SegmentPIDArray,SegmentIDArray,Segment2PIDArray,Segment2IDArray, RecoShowerArrayN)) continue;

                        // Now apply cut conditions: GS  GAMMA SEARCH Alg  --------------------

                        // Check if IP of both to vtx (BT) is ok:
                        IP_Pair_To_InBT_Seg2  =CalcIP(Segment2, vtx);
                        if (IP_Pair_To_InBT_Seg2>CUT_PARAMETER[0]) continue;

                        // if InBT is flagged as MC InBT, take care that only BG or same MC basetracks are taken:
                        if (InBT->MCEvt()>0) if (Segment->MCEvt()>0&&Segment2->MCEvt()>0) if (Segment->MCEvt()!=Segment2->MCEvt()) continue;
                        if (InBT->MCEvt()>0) if (Segment->MCEvt()>0&&Segment2->MCEvt()>0) if (Segment->MCEvt()!=InBT->MCEvt()) continue;
                        if (InBT->MCEvt()>0) if (Segment->MCEvt()>0&&Segment2->MCEvt()>0) if (Segment2->MCEvt()!=InBT->MCEvt()) continue;

                        // In case of two MC events, check for e+ e- pairs
                        // Do this ONLY IF parameter eParaValue[6] is set to choose different Flag() pairs:
                        if (InBT->MCEvt()>0 && CUT_PARAMETER[6]==1) {
                            if (Segment->MCEvt()>0&&Segment2->MCEvt()>0) {
                                if ((Segment2->Flag()+Segment->Flag())!=0) continue;
                            }
                        }

                        // a) Check dR between tracks:
                        if (GetdeltaRWithPropagation(Segment,Segment2)>CUT_PARAMETER[2]) continue;
                        // b) Check dT between tracks:
                        if (GetdeltaThetaSingleAngles(Segment,Segment2)>CUT_PARAMETER[4]) continue;
                        // c) Check dMinDist between tracks:
                        if (GetMinimumDist(Segment,Segment2)>CUT_PARAMETER[1]) continue;

                        // f) Check if this is not a possible fake doublet which is
                        //   sometimes caused by view overlap in the scanning:
                        //    in the EdbPVRQuality class this will be done at start for the whole
                        //    PVR object so this will be later on obsolete.
                        ///if (IsPossibleFakeDoublet(Segment,Segment2) ) continue;
                        //
                        // end of    cut conditions: GS  GAMMA SEARCH Alg  --------------------
                        //


                        if (gEDBDEBUGLEVEL>3) {
                            cout << "EdbShowAlg_GS::FindPairs Pair (PID:" << Segment->PID() << ",ID:" << Segment->ID()<< ";"<< Segment2->PID() << "," << Segment2->ID() << ") has passed all cuts w.r.t to InBT:" << endl;
                            cout << "EdbShowAlg_GS::FindPairs GetdeltaRWithPropagation(Segment,Segment2)  = " << GetdeltaRWithPropagation(Segment,Segment2) << endl;
                            cout << "EdbShowAlg_GS::FindPairs GetdeltaThetaSingleAngles(Segment,Segment2)  = " << GetdeltaThetaSingleAngles(Segment,Segment2) << endl;
                            cout << "EdbShowAlg_GS::FindPairs GetMinimumDist(Segment,Segment2)  = " << GetMinimumDist(Segment,Segment2) << endl;
                            cout << "EdbShowAlg_GS::FindPairs CalcIP(BetterSegment,InBT)  = " << IP_Pair_To_InBT << endl;
                        }

                        if (gEDBDEBUGLEVEL>3)  cout <<"------------"<< endl;

                        // And Add Segment2 to to shower array:
//                         cout << "// Add Segment2 to to shower array (AddBTToArrayWithCeck):" << endl;
                        Bool_t isContained = AddBTToArrayWithCeck(Segment2, GLOBAL_ShowerSegArray);
                        //PrintShowerObjectArray(GLOBAL_ShowerSegArray);
//             cout << " isContained == " << isContained << endl;
                        if (isContained==kTRUE) continue;

                    } //for (Int_t pat_two_bt_cnt=0; ...

                } //for (Int_t pat_one_bt_cnt=0; ...

            } // for (Int_t pat_two_cnt=0; ...


            // Now here do the usual rest for BG density calculation:
            //...
            // Calc BT density around shower:
            EdbPattern* pat_interim=local_gAli->GetPattern(pat_one_cnt);
            CalcTrackDensity(pat_interim,local_gAli_pat_interim_halfsize,npat_int,npat_total,npatN);

            // Calc TrackNumbers for plate for efficency numbers:
            CalcEfficencyNumbers(pat_interim, InBT->MCEvt(), NBT_Neff, NBTMC_Neff,NBTMCe_Neff);

        } //for (Int_t pat_one_cnt=0; ...


        ///============================================================================================
        ///============================================================================================
        ///============================================================================================

        if (gEDBDEBUGLEVEL>2) cout << " GLOBAL_ShowerSegArray->GetEntries() " << GLOBAL_ShowerSegArray->GetEntries() << endl;

        ///----------------------
        /// ONLY FOR THE ANN TRAINING FOR GS ALGO:
        /// DO THIS ONLY IF YOU HAVE PAIRS !!!

        if (GLOBAL_ShowerSegArray->GetEntries()>1) {

            Segment=(EdbSegP*)GLOBAL_ShowerSegArray->At(0);
            Segment2=(EdbSegP*)GLOBAL_ShowerSegArray->At(1);



            if (gEDBDEBUGLEVEL>2) {
                cout << "ONLY FOR THE ANN TRAINING FOR GS ALGO: "  << endl;
                cout << "DO THIS ONLY IF YOU HAVE PAIRS !!!: "  << endl;
                cout << "YES WE HAVE A PAIR Print all segs in the array: "  << endl;
                for (int l=0; l<GLOBAL_ShowerSegArray->GetEntries(); ++l) {
                    EdbSegP* s=(EdbSegP*)GLOBAL_ShowerSegArray->At(l);
                    s->PrintNice();
                }
//      Segment->PrintNice();
//      Segment2->PrintNice();
                cout << "Print again the InBT for crosscheck: InBT= " << endl;
                InBT->PrintNice();
                cout << "Address of Segment   = " << Segment << endl;
                cout << "Address of Segment2  = " << Segment2 << endl;
            }



            // Implicitly, we assume that InBT= Seg->At(0), which I checked is
            // correct. Also in the code it is done this way that InBT is added
            // as first BT.
            //Segment->PrintNice();
            //Segment2->PrintNice();
            //InBT->PrintNice();

            IP_Pair_To_InBT_Seg2=CalcIP(Segment2,vtx);
            IP_Pair_To_InBT=CalcIP(Segment,vtx);

            // Check if both basetracks have a vertex which is upstream
            // of both tracks (only then the two BT are really pointing).
            TObjArray *segments = new TObjArray(2);
            segments->Clear();
            segments->Add(Segment);
            segments->Add(Segment2);
            EdbVertex* vetex = new EdbVertex();
            vetex = CalcVertex(segments);
            cout << "Calculated helper _vetex_ out of the first two segments. " << endl;
            cout <<"vetex ->X,Y,Z:  " << vetex ->X() << " " << vetex ->Y() << " " << vetex ->Z()<< endl;
            if (vetex ->Z()> TMath::Min(Segment->Z(),Segment2->Z()) ) cout << "The interims vertex to which the two segments points lies downstream, or in between the to segmnets, i.e. the do NOT point to an originating vertex. continue now...." << endl;
            cout << "Check:  TMath::Min(Segment->Z(),Segment2->Z())  " <<  TMath::Min(Segment->Z(),Segment2->Z())  << endl;
            cout <<"Segment->Z(): " << Segment->Z() << "  Segment2->Z() " << Segment2->Z() << endl;
            cout << "Address of Segment   = " << Segment << endl;
            cout << "Address of Segment2  = " << Segment2 << endl;
            if (vetex ->Z()> TMath::Min(Segment->Z(),Segment2->Z()) ) {
                cout << " Hmm, do nothing..." << endl;
            }
            else {

                Float_t IP_Seg1ToVtxSeg1Seg2=0;
                Float_t IP_Seg2ToVtxSeg1Seg2=0;
                IP_Seg1ToVtxSeg1Seg2 = CalcIP(Segment,vetex);
                IP_Seg2ToVtxSeg1Seg2 = CalcIP(Segment2,vetex);

                Int_t eRecoMode=2;
                cout << "WARNING   eRecoMode konstant set to =2   TODO ....   change on prompt!!" << endl;
                Float_t eValueGSNN_var00;
                Float_t IP_InBT_To_Vtx=IP_Pair_To_InBT;
                if (eRecoMode==0) eValueGSNN_var00=IP_InBT_To_Vtx;
                if (eRecoMode==1) eValueGSNN_var00=IP_InBT_To_Vtx;
                if (eRecoMode==2) {
                    eValueGSNN_var00=TMath::Min(IP_Seg1ToVtxSeg1Seg2,IP_Seg2ToVtxSeg1Seg2);
                }

                h_GSNN_var00->Fill(eValueGSNN_var00);
                h_GSNN_var01->Fill(GetMinimumDist(Segment,Segment2));
                h_GSNN_var02->Fill(GetdeltaRWithPropagation(Segment,Segment2));
                h_GSNN_var03->Fill(InBT->Z()-vtx->Z());
                h_GSNN_var04->Fill(GetdeltaThetaSingleAngles(Segment,Segment2));
                h_GSNN_var05->Fill(TMath::Abs(Segment->PID()-Segment2->PID()));
                h_GSNN_var06->Fill(Segment2->Flag()+Segment->Flag());

                // Purity 1:    Input =  1.0;
                // Purity 0.5:  Input =  0.5;
                // Purity else: Input =  0.0;
                if (Segment2->Flag()+Segment->Flag()==0&&TMath::Abs(Segment2->Flag())==11&&Segment->MCEvt()>0) {
                    value_GSNN_varInput=1;
                }
                else if (Segment2->Flag()+Segment->Flag()!=0&&TMath::Abs(Segment2->Flag())==11&&Segment->MCEvt()>0) {
                    value_GSNN_varInput=0.5;
                }
                else {
                    value_GSNN_varInput=0;
                }
                value_GSNN_var00=eValueGSNN_var00;
                value_GSNN_var01=GetMinimumDist(Segment,Segment2);
                value_GSNN_var02=GetdeltaRWithPropagation(Segment,Segment2);
                value_GSNN_var03=InBT->Z()-vtx->Z();
                value_GSNN_var04=GetdeltaThetaSingleAngles(Segment,Segment2);
                value_GSNN_var05=TMath::Abs(Segment->PID()-Segment2->PID());
                value_GSNN_var06=Segment2->Flag()+Segment->Flag();

                t_GSNN->Fill();

                cout << "I have filled the GSNN Tree now. End of this bracket." << endl;
            }
        }
        ///----------------------

        if (gEDBDEBUGLEVEL>2) PrintShowerObjectArray(GLOBAL_ShowerSegArray);

        Int_t s_NBT=0;
        Int_t s_NBTMC=0;
        Int_t s_NBTallMC=0;
        Int_t s_NBTeMC=0;
        Double_t  s_eff=0;
        Double_t s_purall=0;
        Double_t s_pure=0;
        CalcEffPurOfShower(GLOBAL_ShowerSegArray, s_NBT, s_NBTMC, s_NBTallMC, s_NBTeMC, s_purall, s_pure);

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
//         cout << "// 4) Calculate pur/eff/NBT numbers...done." << endl;


//     XXXDEBUG

        //-----------------------------------
        // 5) Fill Tree:
        //-----------------------------------
        if (gEDBDEBUGLEVEL>2) cout << "// 5) Fill Tree:" << endl;
        TREE_ShowRecEff->Fill();
        if (gEDBDEBUGLEVEL>3) TREE_ShowRecEff->Show(TREE_ShowRecEff->GetEntries()-1);


        //-----------------------------------
        // 6a) Transfer ShowerArray to treebranchTreeEntry:
        //-----------------------------------
        if (cmd_OUTPUTLEVEL>0) {
            //cout << "// 6a) Transfer ShowerArray to treebranchTreeEntry:" << endl;
            TransferShowerObjectArrayIntoEntryOfTreebranchShowerTree(TREE_ShowShower,GLOBAL_ShowerSegArray);
        }


        //------------------------------------
        // Reset and delete important things:
        // also  to avoid memory problems ...
        //-----------------------------------
        GLOBAL_ShowerSegArray->Clear();
        if (gEDBDEBUGLEVEL>2) cout << "--- --- after clear:  GLOBAL_ShowerSegArray->GetEntries(): "<< GLOBAL_ShowerSegArray->GetEntries() << endl;
        delete local_gAli;
        local_gAli=0;
        ++GLOBAL_INBTSHOWERNR;
        //------------------------------------



        IsFirstLoopCount=kFALSE;
    }
    // end of loop over GLOBAL_InBTArrayEntries
    //-----------------------------------------------------------------

    if (gEDBDEBUGLEVEL==2) cout << endl<<flush;
    if (gEDBDEBUGLEVEL>3) cout << "---TREE_ShowRecEff->GetEntries() ... " << TREE_ShowRecEff->GetEntries() << endl;
    if (gEDBDEBUGLEVEL>3) cout << "---GLOBAL_INBTSHOWERNR ... " << GLOBAL_INBTSHOWERNR<< endl;



    Log(2, "ShowRec.cpp", "--- void ReconstructShowers_GS() done.");
    return;
}
//-------------------------------------------------------------------------------------------
