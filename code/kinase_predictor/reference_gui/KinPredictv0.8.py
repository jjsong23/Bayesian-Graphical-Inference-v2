# -*- coding: utf-8 -*-
"""
Created on Thu May  5 10:33:44 2022

@author: Brian
"""
# AQP2 Ser256: EVRRRQSVELHSP

import sys
import glob
import os
from operator import itemgetter
from pathlib import Path
import numpy as py
import pandas as pd
from KinPredUITable4 import Ui_Form
from PyQt5.QtWidgets import QApplication, QWidget, QFileDialog, QMessageBox, QTableWidgetItem

class MyApp(QWidget):
    
    def __init__(self):
        super().__init__()
        self.ui = Ui_Form()
        self.ui.setupUi(self)
        self.setWindowTitle('KinasePredictor')
        self.ui.seqsubmit.clicked.connect(self.kinpred)
        self.ui.clearinput.clicked.connect(self.ui.userinput.clear)
        self.ui.clearoutput.clicked.connect(self.tableclear)
        self.ui.browse.clicked.connect(self.browsefiles)
        self.ui.exportoutput.clicked.connect(self.exportresults)
        self.ui.kinase_probability=[]
        self.ui.exportfile=[]
        

    def kinpred(self):
        userseq=self.ui.userinput.text()
        rowCount=self.ui.outputtable.rowCount()
        if len(userseq)!=13:
            # def show_popup(self):
                msg=QMessageBox()
                msg.setWindowTitle("ERROR")
                msg.setText("Input sequence must be 13 amino acids long")
                x = msg.exec_()
                return
            # self.ui.outputtable.insertRow(rowCount)
        else:
            self.usersequpper=userseq.upper()
        # filepath to kinase output matrices
        path=Path('./output_matrices/Ser-Thr_output_matrices/' )
        path2=Path('./output_matrices/Tyr_output_matrices/' )
        # create list of all kinase output matrices in filepath
        if self.ui.radioButton.isChecked():
            csv_files = glob.glob(os.path.join(path2, "*.csv"))
        elif self.ui.radioButton_2.isChecked():
            csv_files = glob.glob(os.path.join(path, "*.csv"))
        else:
            msg=QMessageBox()
            msg.setWindowTitle("ERROR")
            msg.setText("Select a Center Amino Acid (Y or S/T)")
            x = msg.exec_()
            
        # csv_files = glob.glob(os.path.join(path, "*.csv"))
        # csv_files2 = glob.glob(os.path.join(path2, "*.csv"))
        # define row number for each amino acid in kinase output matrix files
        aarow={'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 'K': 8, 'L': 9, 'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14,
                'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19, 'B': 20, 'J':21}
        # create lists for coordinates of input sequence and information content for each kinase at those coords
        kincoord=[]
        kininfolist=[]
        # kinproblist=[]
        # generate coordinate list for user input sequence, correspond to AA and position in output matrices
        for pos,AA in enumerate(self.usersequpper):
            x=aarow[AA]
            y=pos+1
            kincoord+=[[x,y]]

        # Iterate through each kinase csv file, get the value at each coordinate, and add these values for scalar kinase product
            
        # for i,j  in zip(csv_files,csv_files2):
        for i  in csv_files:
            df=pd.read_csv(i)
            # dg=pd.read_csv(j)
        # Pair each value with the corresponding kinase filename
            kinase=os.path.basename(i).split('.')[0]
            # kinase2=os.path.basename(j).split('.')[0]
            kinprob=[]
            # kinprob2=[]
            for k,l in kincoord:
                aainfo=df.iloc[k][l]
                # aaprob=dg.iloc[k][l]
                kinprob+=[aainfo]
                # kinprob2+=[aaprob]
        # Append each kinase value pair to a list, sort that list in descending order, and display top 10 kinase values.
            kinpred=round(py.sum(kinprob),4)
            # kinpred2=round(py.sum(kinprob2),4)
            kininfolist+=[[kinase,kinpred]]
            # kinproblist+=[[kinase2,kinpred2]]
            # self.ui.outputtable.insertRow(rowCount)
            # self.ui.outputtable.setItem(rowCount-2, 4, QTableWidgetItem(kinase))
            # self.ui.outputtable.setItem(rowCount-2, 5, QTableWidgetItem(str(kinpred)))
            # self.ui.outputtable.setItem(rowCount-1, 6, QTableWidgetItem(str(kinase2)))
            # self.ui.outputtable.setItem(rowCount-1, 7, QTableWidgetItem(str(kinpred2)))
        self.ui.outputtable.setSortingEnabled(True)
        kinoutputsort=pd.DataFrame(sorted(kininfolist, key=itemgetter(1), reverse=False), columns= ['Kinase', 'Information Content'])
        kininfosort=pd.DataFrame(sorted(kininfolist, key=itemgetter(1), reverse=True), columns= ['Kinase', 'Information Content'])
        for i in range(len(kinoutputsort)):
            self.ui.outputtable.insertRow(rowCount)
            self.ui.outputtable.setItem(rowCount-2, 4, QTableWidgetItem(kinoutputsort.iloc[i]['Kinase']))
            self.ui.outputtable.setItem(rowCount-2, 5, QTableWidgetItem(str(kinoutputsort.iloc[i]['Information Content'])))
            # self.ui.outputtable.setItem(rowCount-2, 4, QTableWidgetItem(kininfosort.iloc[i]['Kinase']))
            # self.ui.outputtable.setItem(rowCount-2, 5, QTableWidgetItem(str(kininfosort.iloc[i]['Information Content'])))
        # kinprobsort=pd.DataFrame(sorted(kinproblist, key=itemgetter(1), reverse=True), columns= ['Kinase2','Probability Sum'])
        # self.kinase_probability=pd.concat([kininfosort, kinprobsort], axis=1)
        self.kinase_probability=pd.concat([kininfosort], axis=1)
        # listToStr = '\n'.join(map(str, self.kinase_probability))
        # self.ui.outputwindow.setText(listToStr)
        return self.usersequpper
        return self.kinase_probability
        return kincoord
    
  #  exportfile = kinpred(self)
    
    def browsefiles(self):
        fpath=QFileDialog.getSaveFileName(self,'Open File', 'C:\\Users\\*\\Desktop\\', 'CSV files (*.csv)')
        self.ui.filename.setText(fpath[0])
    
    def exportresults(self):
        fname=self.ui.filename.text()
        if fname=='':
            msg2=QMessageBox()
            msg2.setWindowTitle("ERROR")
            msg2.setText("Designate a folder and filename before exporting")
            x = msg2.exec_()
        else:
            df = self.kinase_probability
            # df = pd.DataFrame(self.kinase_probability, columns= ['Kinase', 'Information Content'])
            df.index.name = "Input Sequence: "+self.usersequpper
            df.to_csv(fname)
            msg=QMessageBox()
            msg.setWindowTitle("Export Results")
            msg.setText("Export Complete")
            x = msg.exec_()
    # def exportresults(self):
    #     listToStr=kinpred(self)
    #     listToStr.to_csv('output.csv') 
    def tableclear(self):
        self.ui.outputtable.clearContents
        self.ui.outputtable.setRowCount(0)
        self.ui.outputtable.setSortingEnabled(False)
                
if __name__ == '__main__':
    app=QApplication(sys.argv)
    myApp=MyApp()
    myApp.show()
    
    try:
        sys.exit(app.exec())
    except SystemExit:
        print('Closing Window...')
        

# import glob
# import os
# from operator import itemgetter
# import numpy as py
# import pandas as pd

  #  self.browse.clicked.connect(self.browsefiles)
    
  #  def browsefiles(self):
  #      fname=QFileDialog.getOpenFileName(self,'Open File', 'C:\Users\*\Desktop', 'CSV files (*.csv)')
  #      self.filename.setText(fname[0])
        
  # path='H:\Kinase Logo Website\Kinase Prediction tool\output_matrices\output_matrices'

# def kinpred(userseq): 
#     # filepath to kinase output matrices
#     path='H:\Kinase Logo Website\Kinase Prediction tool\output_matrices\output_matrices'
#     # create list of all kinase output matrices in filepath
#     csv_files = glob.glob(os.path.join(path, "*.csv"))
#     # define row number for each amino acid in kinase output matrix files
#     aarow={'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 'K': 8, 'L': 9, 'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14,
#        'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19, 'B': 20, 'J':21}
#     # create lists for coordinates of input sequence and information content for each kinase at those coords
#     kincoord=[]
#     kinproblist=[]
#     # generate coordinate list for user input sequence, correspond to AA and position in output matrices
#     for pos,AA in enumerate(userseq):
#         x=aarow[AA]
#         y=pos+1
#         kincoord+=[[x,y]]
#     print(kincoord)
#     # Iterate through each kinase csv file, get the value at each coordinate, and add these values for scalar kinase product
#     for i in csv_files:
#         df=pd.read_csv(i)
#         # Pair each value with the corresponding kinase filename
#         kinase=os.path.basename(i).split('.')[0]
#         kinprob=[]
#         for i,j in kincoord:
#             aainfo=df.iloc[i][j]
#             kinprob+=[aainfo]
#     # Append each kinase value pair to a list, sort that list in descending order, and display top 10 kinase values.
#         kinpred=py.sum(kinprob)
#         kinproblist+=[[kinase,kinpred]]
#     kinase_probability=sorted(kinproblist, key=itemgetter(1), reverse=True)
#     toptenkinase=kinase_probability[0:10]
#     print(toptenkinase)
