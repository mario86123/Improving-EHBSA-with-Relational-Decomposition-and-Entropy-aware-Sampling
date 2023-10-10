CC = g++ -O2 -Wno-deprecated -std=c++11

all: RankingEDAsCEC

RankingEDAsCEC: main.o PBP.o PFSP.o FSSP.o LOP.o QAP.o TSP.o SOP.o API.o SDF.o GAP.o VRP.o Tools.o Cayley.o Distance.o Kendall.o GeneralizedKendall.o GeneralizedDistance.o RankingModel.o GeneralizedMallowsModel.o MallowsModel.o BWModel.o MRModel.o EHBSAWO.o NHBSAWO.o EHBSAWT.o ENEHBSAWT.o NHBSAWT.o TEMP.o RRBOP.o FERBOP.o ARBOP.o PRBOP.o NRBOP.o ANRBOP.o ACNRBOP.o CHI.o ACHIGN.o AGNRBOP.o PREPRUNE.o FOUR.o SRBOP.o FCRBOP.o RBOP.o BM.o NNHBSAWT.o EEHBSAWT.o SNHBSAWT.o  MRWT.o AMRModel.o FMSTModel.o LLMST.o MST2.o MST2ME.o MST3.o MST4.o SoftMST.o RMST.o SRMST.o AMST.o MSTModel.o MSTME.o Population.o Individual.o RankingEDA.o GeneralizedCayley.o Ferrers_diagram.o Ferrers_diagram2.o Generic.o Ulam.o Ulam2.o NewtonRaphson.o
	$(CC)  main.o PBP.o PFSP.o FSSP.o LOP.o QAP.o TSP.o SOP.o API.o SDF.o GAP.o VRP.o Tools.o Cayley.o Distance.o Kendall.o GeneralizedKendall.o GeneralizedDistance.o RankingModel.o GeneralizedMallowsModel.o MallowsModel.o BWModel.o MRModel.o EHBSAWO.o NHBSAWO.o EHBSAWT.o ENEHBSAWT.o NHBSAWT.o TEMP.o RRBOP.o FERBOP.o ARBOP.o PRBOP.o NRBOP.o ANRBOP.o ACNRBOP.o CHI.o ACHIGN.o AGNRBOP.o PREPRUNE.o FOUR.o SRBOP.o FCRBOP.o RBOP.o BM.o NNHBSAWT.o EEHBSAWT.o SNHBSAWT.o  MRWT.o AMRModel.o FMSTModel.o LLMST.o MST2.o MST2ME.o MST3.o MST4.o SoftMST.o RMST.o SRMST.o AMST.o MSTModel.o MSTME.o Population.o Individual.o RankingEDA.o GeneralizedCayley.o Ferrers_diagram.o Ferrers_diagram2.o Generic.o Ulam.o Ulam2.o NewtonRaphson.o -lstdc++ -lm -o RankingEDAsCEC

main.o: main.cpp RankingEDA.h PBP.h PFSP.h FSSP.h LOP.h QAP.h TSP.h SOP.h API.h SDF.h GAP.h VRP.h
	$(CC) -c main.cpp
 
Cayley.o: Cayley.cpp Cayley.h Tools.h Distance.h Variables.h NewtonRaphson.h
	$(CC) -c Cayley.cpp

Kendall.o: Kendall.cpp Kendall.h Tools.h Distance.h Variables.h NewtonRaphson.h
	$(CC) -c Kendall.cpp
	
Ulam.o: Ulam.cpp Ulam.h Tools.h Distance.h Variables.h Ferrers_diagram.h NewtonRaphson.h
	$(CC) -c Ulam.cpp
	
Ulam2.o: Ulam2.cpp Ulam2.h Tools.h Distance.h Variables.h Ferrers_diagram2.h NewtonRaphson.h
	$(CC) -c Ulam2.cpp
	
Ferrers_diagram.o: Ferrers_diagram.cpp Ferrers_diagram.h Generic.h
	$(CC) -c Ferrers_diagram.cpp

Ferrers_diagram2.o: Ferrers_diagram2.cpp Ferrers_diagram2.h Generic.h
	$(CC) -c Ferrers_diagram2.cpp
	
Generic.o: Generic.cpp Generic.h
	$(CC) -c Generic.cpp

Distance.o: Distance.cpp Distance.h 
	$(CC) -c Distance.cpp

GeneralizedKendall.o: GeneralizedKendall.cpp GeneralizedKendall.h Tools.h GeneralizedDistance.h NewtonRaphson.h
	$(CC) -c GeneralizedKendall.cpp

GeneralizedCayley.o: GeneralizedCayley.cpp GeneralizedCayley.h Tools.h GeneralizedDistance.h NewtonRaphson.h
	$(CC) -c GeneralizedCayley.cpp

GeneralizedDistance.o: GeneralizedDistance.cpp GeneralizedDistance.h
	$(CC) -c GeneralizedDistance.cpp

RankingModel.o: RankingModel.cpp RankingModel.h Population.h
	$(CC) -c RankingModel.cpp

MallowsModel.o: MallowsModel.cpp MallowsModel.h Individual.h RankingModel.h Population.h Distance.h Kendall.h Cayley.h Tools.h Ulam.h Ulam2.h
	$(CC) -c MallowsModel.cpp

GeneralizedMallowsModel.o: GeneralizedMallowsModel.cpp GeneralizedMallowsModel.h Individual.h RankingModel.h Population.h GeneralizedCayley.h  GeneralizedDistance.h GeneralizedKendall.h Tools.h
	$(CC) -c GeneralizedMallowsModel.cpp

BWModel.o: BWModel.cpp BWModel.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c BWModel.cpp

MRModel.o: MRModel.cpp MRModel.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c MRModel.cpp

EHBSAWO.o: EHBSAWO.cpp EHBSAWO.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c EHBSAWO.cpp

EHBSAWT.o: EHBSAWT.cpp EHBSAWT.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c EHBSAWT.cpp

ENEHBSAWT.o: ENEHBSAWT.cpp ENEHBSAWT.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c ENEHBSAWT.cpp

NHBSAWO.o: NHBSAWO.cpp NHBSAWO.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c NHBSAWO.cpp

NHBSAWT.o: NHBSAWT.cpp NHBSAWT.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c NHBSAWT.cpp

TEMP.o: TEMP.cpp TEMP.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c TEMP.cpp

RRBOP.o: RRBOP.cpp RRBOP.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c RRBOP.cpp

SRBOP.o: SRBOP.cpp SRBOP.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c SRBOP.cpp

RBOP.o: RBOP.cpp RBOP.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c RBOP.cpp

BM.o: BM.cpp BM.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c BM.cpp

ARBOP.o: ARBOP.cpp ARBOP.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c ARBOP.cpp

PRBOP.o: PRBOP.cpp PRBOP.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c PRBOP.cpp

NRBOP.o: NRBOP.cpp NRBOP.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c NRBOP.cpp

ANRBOP.o: ANRBOP.cpp ANRBOP.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c ANRBOP.cpp

ACNRBOP.o: ACNRBOP.cpp ACNRBOP.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c ACNRBOP.cpp

CHI.o: CHI.cpp CHI.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c CHI.cpp

ACHIGN.o: ACHIGN.cpp ACHIGN.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c ACHIGN.cpp

AGNRBOP.o: AGNRBOP.cpp AGNRBOP.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c AGNRBOP.cpp
	
PREPRUNE.o: PREPRUNE.cpp PREPRUNE.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c PREPRUNE.cpp

FOUR.o: FOUR.cpp FOUR.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c FOUR.cpp

FERBOP.o: FERBOP.cpp FERBOP.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c FERBOP.cpp

FCRBOP.o: FCRBOP.cpp FCRBOP.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c FCRBOP.cpp

NNHBSAWT.o: NNHBSAWT.cpp NNHBSAWT.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c NNHBSAWT.cpp

EEHBSAWT.o: EEHBSAWT.cpp EEHBSAWT.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c EEHBSAWT.cpp

SNHBSAWT.o: SNHBSAWT.cpp SNHBSAWT.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c SNHBSAWT.cpp

MRWT.o: MRWT.cpp MRWT.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c MRWT.cpp

AMRModel.o: AMRModel.cpp AMRModel.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c AMRModel.cpp

FMSTModel.o: FMSTModel.cpp FMSTModel.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c FMSTModel.cpp

MST2.o: MST2.cpp MST2.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c MST2.cpp

MST2ME.o: MST2ME.cpp MST2ME.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c MST2ME.cpp

MST3.o: MST3.cpp MST3.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c MST3.cpp

MST4.o: MST4.cpp MST4.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c MST4.cpp

SoftMST.o: SoftMST.cpp SoftMST.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c SoftMST.cpp

LLMST.o: LLMST.cpp LLMST.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c LLMST.cpp

RMST.o: RMST.cpp RMST.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c RMST.cpp

SRMST.o: SRMST.cpp SRMST.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c SRMST.cpp

AMST.o: AMST.cpp AMST.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c AMST.cpp

MSTModel.o: MSTModel.cpp MSTModel.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c MSTModel.cpp

MSTME.o: MSTME.cpp MSTME.h Individual.h RankingModel.h Population.h Tools.h
	$(CC) -c MSTME.cpp

Population.o: Population.cc Population.h Individual.h Tools.h
	$(CC) -c Population.cc

Individual.o: Individual.cc Individual.h Variables.h
	$(CC) -c Individual.cc

RankingEDA.o: RankingEDA.cpp RankingEDA.h MallowsModel.h GeneralizedMallowsModel.h BWModel.h Tools.h PBP.h Population.h RankingModel.h
	$(CC) -c RankingEDA.cpp

Tools.o: Tools.cpp Tools.h
	$(CC) -c Tools.cpp

PBP.o: PBP.cpp PBP.h
	$(CC) -c PBP.cpp

PFSP.o: PFSP.cpp PFSP.h Tools.h PBP.h
	$(CC) -c PFSP.cpp

FSSP.o: FSSP.cpp FSSP.h Tools.h PBP.h
	$(CC) -c FSSP.cpp

LOP.o: LOP.cpp LOP.h Tools.h PBP.h                        
	$(CC) -c LOP.cpp

QAP.o: QAP.cpp QAP.h Tools.h PBP.h                        
	$(CC) -c QAP.cpp

TSP.o: TSP.cpp TSP.h Tools.h PBP.h                        
	$(CC) -c TSP.cpp

SOP.o: SOP.cpp SOP.h Tools.h PBP.h                        
	$(CC) -c SOP.cpp

API.o: API.cpp API.h Tools.h PBP.h                        
	$(CC) -c API.cpp

SDF.o: SDF.cpp SDF.h Tools.h PBP.h                        
	$(CC) -c SDF.cpp

GAP.o: GAP.cpp GAP.h Tools.h PBP.h                        
	$(CC) -c GAP.cpp

VRP.o: VRP.cpp VRP.h Tools.h PBP.h                        
	$(CC) -c VRP.cpp
	
NewtonRaphson.o: NewtonRaphson.cpp NewtonRaphson.h                         
	$(CC) -c NewtonRaphson.cpp

clean:
	rm -f core *.o RankingEDAsCEC

