# -*- coding: utf-8 -*-

import torch
from torch import nn
import torch.nn.functional as F
import numpy as np
import random
from sklearn.metrics import roc_auc_score
from sklearn import metrics
from sklearn.metrics import accuracy_score
from sklearn.metrics import f1_score
import os
import datetime

import src.config as config
from src.load_dataset import load_cross_data,load_test_data,load_train_data
from src.load_dataset import load_cross_data_for_fivefold,load_test_data_for_fivefold,load_train_data_for_fivefold
import time

def train_model_for_human(dataset_name,split_mode,ontology_name,iea_name,part_name):
    train_start_time = time.time()
    dataset_name=dataset_name
    split_mode=split_mode #"connectivity"
    #ontology_name='BiologicalProcess'
    ontology_name=ontology_name #'CellularComponent'
    #ontology_name='MolecularFunction'
    #torch.manual_seed(1)    # reproducible
    EPOCH = 5

    # Hyper Parameters
      # train the training data n times, to save time, we just train 1 epoch
    BATCH_SIZE = 256  # rnn time step
    GOdim = config.GODIM
    INPUT_SIZE = GOdim  # rnn input size
    LR = 0.01  # learning rate

    #
    #
    # if config.IEA_FLAG == True:
    #     iea_name = "iea+"
    # else:
    #     iea_name = "iea-"

    def createFolder(directory):
        try:
            if not os.path.exists(directory):
                os.makedirs(directory)
        except OSError:
            print('Error: Creating directory. ' + directory)


    nowtime = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")

    if split_mode=="fivefold":
        LogFileFolder1 = config.home_path + '/train_log/fivefold/' +dataset_name+"/"+ 'GOdim' + str(GOdim) + '_' + iea_name + '_' + ontology_name
        LogFileFolder2 = LogFileFolder1 + '/'+nowtime+part_name
    else:
        LogFileFolder1 = config.home_path + '/train_log/connectivity/'+dataset_name + 'GOdim' + str(
            GOdim) + '_' + iea_name + '_' + ontology_name
        LogFileFolder2 = LogFileFolder1 + '/' + nowtime

    createFolder(LogFileFolder1)
    createFolder(LogFileFolder2)

    AUCList = []
    AccuracyList = []
    AUCList.append(-1)
    AccuracyList.append(-1)

    path = config.home_path
    model_path = os.path.join(LogFileFolder2, 'model')
    pred_cross_score_file = os.path.join(LogFileFolder2, 'CrossPreditionScore')
    pred_test_score_file=os.path.join(LogFileFolder2, 'TestPreditionScore')
    Cross_Result_Log_File = open(os.path.join(LogFileFolder2, 'CrossResultListLog'), 'w')
    Test_Result_Log_File = open(os.path.join(LogFileFolder2, 'TestResultListLog'), 'w')
    Train_Log_File = open(os.path.join(LogFileFolder2, 'TrainLog'), 'w')

    if split_mode=="connectivity":
        # 加载训练集
        p1_train, p2_train, label_train =load_train_data(dataset_name=dataset_name,ontology_name=ontology_name,split_mode=split_mode,iea_name=iea_name)
        # 加载验证集
        p1_cross, p2_cross, label_cross =load_cross_data(dataset_name=dataset_name,ontology_name=ontology_name,split_mode=split_mode,iea_name=iea_name)

        # 加载测试集
        p1_test, p2_test, label_test = load_test_data(dataset_name=dataset_name, ontology_name=ontology_name,
                                                          split_mode=split_mode, iea_name=iea_name)
    else:
        # 加载训练集
        p1_train, p2_train, label_train = load_train_data_for_fivefold(dataset_name=dataset_name, ontology_name=ontology_name, split_mode=split_mode,
                                        iea_name=iea_name,part_name=part_name)
        # 加载验证集
        p1_cross, p2_cross, label_cross = load_cross_data_for_fivefold(dataset_name=dataset_name, ontology_name=ontology_name,
                                                       split_mode=split_mode, iea_name=iea_name,part_name=part_name)
        # 加载测试集
        p1_test, p2_test, label_test = load_test_data_for_fivefold(dataset_name=dataset_name,
                                                                       ontology_name=ontology_name,
                                                                       split_mode=split_mode, iea_name=iea_name,
                                                                       part_name=part_name)

    len_train=len(p1_train)
    print('train data length {}'.format(len_train))
    len_cross = len(p1_cross)
    print('cross data length {}'.format(len_cross))
    len_test = len(p1_test)
    print('test data length {}'.format(len_test))

    class RNN(nn.Module):
        def __init__(self):
            super(RNN, self).__init__()

            self.rnn = nn.LSTM(  # if use nn.RNN(), it hardly learns
                input_size=INPUT_SIZE,
                hidden_size=64,  # 蛋白质向量维度
                num_layers=1,  # number of rnn layer
                batch_first=True,  # input & output will has batch size as 1s dimension. e.g. (batch, time_step, input_size)
            )

            # self.out = nn.Linear(64, 10)
            self.logic_hidden1 = torch.nn.Linear(128, 64)  # （蛋白质向量维度*2，
            self.logic_hidden2 = torch.nn.Linear(64, 16)  # 隐藏层线性输出
            self.logic_hidden3 = torch.nn.Linear(16, 8)  # 隐藏层线性输出
            self.logic_out = torch.nn.Linear(8, 2)  # 输出层线性输出

        def forward(self, x1, x2):
            # x shape (batch, time_step, input_size)
            # r_out shape (batch, time_step, output_size)
            # h_n shape (n_layers, batch, hidden_size)
            # h_c shape (n_layers, batch, hidden_size)
            r_out1, (h_n1, h_c1) = self.rnn(x1, None)  # None represents zero initial hidden state

            # choose r_out at the last time step
            out1 = r_out1[:, -1, :]
            r_out2, (h_n2, h_c2) = self.rnn(x2, None)  # None represents zero initial hidden state

            # choose r_out at the last time step
            out2 = r_out2[:, -1, :]
            x = torch.cat((out1, out2), 1)
            # x=torch.transpose(x, 0, 1)
            x = F.relu(self.logic_hidden1(x))
            x = F.relu(self.logic_hidden2(x))
            x = F.relu(self.logic_hidden3(x))
            x = self.logic_out(x)
            return x, out1, out2


    # new_rnn=RNN()
    # new_rnn.load_state_dict(torch.load(model_path))
    # trian(new_rnn)


    def train(rnn):
        optimizer = torch.optim.Adam(rnn.parameters(), lr=LR)  # optimize all cnn parameters
        loss_func = nn.CrossEntropyLoss()
        n = len(p1_train)
        best_cross_model=None #交叉验证获得的最好的模型
        # print('Starting training................')
        for epoch in range(EPOCH):
            print('epoch')
            print(epoch)
            Train_Log_File.write("epoch:  " + str(epoch) + '\n')
            # random.shuffle(p1)
            for st in range(0, n, BATCH_SIZE):
                # create batch data
                en = st + BATCH_SIZE
                if st > n - 1:
                    break
                if en > n:
                    en = n
                print("st:  " + str(st) + "en :" + str(en))

                loss = 0
                for step in range(st, en):  # batch train

                    x1 = torch.from_numpy(p1_train[step][0]).type(torch.FloatTensor)
                    x2 = torch.from_numpy(p2_train[step][0]).type(torch.FloatTensor)
                    t = torch.from_numpy(np.array([label_train[step]])).type(torch.LongTensor)
                    for i in range(1, len(p1_train[step])):
                        x1_tmp = torch.from_numpy(p1_train[step][i]).type(torch.FloatTensor)
                        x1 = torch.cat((x1, x1_tmp), 0).type(torch.FloatTensor)
                    for j in range(1, len(p2_train[step])):
                        x2_tmp = torch.from_numpy(p2_train[step][j]).type(torch.FloatTensor)
                        x2 = torch.cat((x2, x2_tmp), 0).type(torch.FloatTensor)
                    x1 = x1.view(1, len(p1_train[step]), GOdim).type(torch.FloatTensor)
                    x2 = x2.view(1, len(p2_train[step]), GOdim).type(torch.FloatTensor)
                    output, p1_train_vector, p2_train_vector = rnn(x1, x2)

                    # simValue=F.cosine_similarity(p1_train_vector,p2_vector).data.numpy().squeeze()
                    # print(simValue)
                    # pred_score=F.softmax(output).data.numpy().squeeze()[1]
                    # print(pred_score)
                    # print(t)
                    loss_tmp = loss_func(output, t)
                    loss = loss + loss_tmp
                    # labeltmp=torch.randint(0,1, (1,), dtype=torch.long)
                    # print('output-real'+str(output)+'--'+str(t))
                print('Average loss')
                print(loss / BATCH_SIZE)
                Train_Log_File.write("Step " + str(step) + ":  " + str(loss / BATCH_SIZE))
                optimizer.zero_grad()  # clear gradients for this training

                loss.backward()
                optimizer.step()
            if epoch % 1 == 0:
                prediction_scores, label_test, auc, accuracy,best_cross_model = model_cross(rnn, loss_func, epoch,best_cross_model)

            # 用best_cross_model 测试数据得到结果
        model_test(best_cross_model, loss_func)

    def model_test(rnn, loss_func):
        test_start_time = time.time()
        TEST_SIZE = len(p1_test)
        loss = 0
        sump = 0
        prediction_scores = []  # 概率
        prediction_value = []  # 真实值
        prediction_similarity = []
        for step in range(TEST_SIZE):  # batch train
            x1 = torch.from_numpy(p1_test[step][0]).type(torch.FloatTensor)
            x2 = torch.from_numpy(p2_test[step][0]).type(torch.FloatTensor)
            t = torch.from_numpy(np.array([label_test[step]])).type(torch.LongTensor)
            for i in range(1, len(p1_test[step])):
                x1_tmp = torch.from_numpy(p1_test[step][i]).type(torch.FloatTensor)
                x1 = torch.cat((x1, x1_tmp), 0).type(torch.FloatTensor)
            for j in range(1, len(p2_test[step])):
                x2_tmp = torch.from_numpy(p2_test[step][j]).type(torch.FloatTensor)
                x2 = torch.cat((x2, x2_tmp), 0).type(torch.FloatTensor)
            x1 = x1.view(1, len(p1_test[step]), GOdim).type(torch.FloatTensor)
            x2 = x2.view(1, len(p2_test[step]), GOdim).type(torch.FloatTensor)

            output, p1_vector, p2_vector = rnn(x1, x2)
            simValue = F.cosine_similarity(p1_vector, p2_vector).data.numpy().squeeze()
            prediction_similarity.append(simValue)
            loss_tmp = loss_func(output, t)
            loss = loss + loss_tmp
            # print('test output realLabel'+str(output)+'--'+str(t))

            # 预测为正样本的概率
            pred_score = F.softmax(output).data.numpy().squeeze()[1]
            # print('pred_score')
            # print(pred_score)
            # print(t)
            prediction_scores.append(pred_score)

            # 分类
            prediction = torch.max(F.softmax(output), 1)[1]
            pred_value = prediction.data.numpy().squeeze()
            # print('pred_value')
            # print(pred_value)
            prediction_value.append(pred_value)
        test_end_time = time.time()
        print('Average test loss')
        ave_losss_test = loss / TEST_SIZE
        print(ave_losss_test)

        print('Accuracy')
        accuracy = accuracy_score(label_test, prediction_value)
        print(accuracy)

        print('F1 socre')
        F1Score = f1_score(label_test, prediction_value, average="macro")
        print(F1Score)

        print('AUC')
        auc = AUC(label_test, prediction_scores)
        print(auc)

        auc_sim = AUC(label_test, prediction_similarity)
        print(auc_sim)
        Test_Result_Log_File.write('Average test loss: ')
        Test_Result_Log_File.write(str(ave_losss_test) + '\n')

        Test_Result_Log_File.write('AUC: ')
        Test_Result_Log_File.write(str(auc) + '\n')

        Test_Result_Log_File.write('AUC_SIM: ')
        Test_Result_Log_File.write(str(auc_sim) + '\n')

        Test_Result_Log_File.write('Accuracy: ')
        Test_Result_Log_File.write(str(accuracy) + '\n')

        Test_Result_Log_File.write('F1: ')
        Test_Result_Log_File.write(str(F1Score) + '\n')

        test_used_time = test_end_time - test_start_time
        Test_Result_Log_File.write("\n")
        Test_Result_Log_File.write("Test Start Time:" + str(time.asctime(time.localtime(test_start_time)))+"\n")
        Test_Result_Log_File.write("Test End Time:" + str(time.asctime(time.localtime(test_end_time)))+"\n")
        Test_Result_Log_File.write("Test Used Time:" + str(test_used_time)+"\n")

        Test_Result_Log_File.flush()
        f = open(pred_test_score_file, 'w')
        for i in range(len(prediction_scores)):
            f.write(str(prediction_scores[i]) + ' ' + str(prediction_similarity[i]) + ' ' + str(label_test[i]) + '\n')
        f.close()
        return prediction_scores, label_test, auc, accuracy

    def model_cross(rnn, loss_func, epoch,best_cross_model):
        CROSS_SIZE = len(p1_cross)
        loss = 0
        sump = 0
        prediction_scores = []  # 概率
        prediction_value = []  # 真实值
        presiction_similarity = []
        for step in range(CROSS_SIZE):  # batch train
            x1 = torch.from_numpy(p1_cross[step][0]).type(torch.FloatTensor)
            x2 = torch.from_numpy(p2_cross[step][0]).type(torch.FloatTensor)
            t = torch.from_numpy(np.array([label_cross[step]])).type(torch.LongTensor)
            for i in range(1, len(p1_cross[step])):
                x1_tmp = torch.from_numpy(p1_cross[step][i]).type(torch.FloatTensor)
                x1 = torch.cat((x1, x1_tmp), 0).type(torch.FloatTensor)
            for j in range(1, len(p2_cross[step])):
                x2_tmp = torch.from_numpy(p2_cross[step][j]).type(torch.FloatTensor)
                x2 = torch.cat((x2, x2_tmp), 0).type(torch.FloatTensor)
            x1 = x1.view(1, len(p1_cross[step]), GOdim).type(torch.FloatTensor)
            x2 = x2.view(1, len(p2_cross[step]), GOdim).type(torch.FloatTensor)

            output, p1_vector, p2_vector = rnn(x1, x2)
            simValue = F.cosine_similarity(p1_vector, p2_vector).data.numpy().squeeze()
            presiction_similarity.append(simValue)
            loss_tmp = loss_func(output, t)
            loss = loss + loss_tmp
            # print('test output realLabel'+str(output)+'--'+str(t))

            # 预测为正样本的概率
            pred_score = F.softmax(output).cpu().data.numpy().squeeze()[1]
            prediction_scores.append(pred_score)
            # 分类
            prediction = torch.max(F.softmax(output), 1)[1]
            pred_value = prediction.cpu().data.numpy().squeeze()
            # print('pred_value')
            # print(pred_value)
            prediction_value.append(pred_value)

        print('Average cross loss')
        ave_losss_cross = loss / CROSS_SIZE
        print(ave_losss_cross)

        print('Accuracy')
        accuracy = accuracy_score(label_cross, prediction_value)
        print(accuracy)

        print('F1 socre')
        F1Score = f1_score(label_cross, prediction_value, average="macro")
        print(F1Score)

        print('AUC')
        auc = AUC(label_cross, prediction_scores)
        print(auc)
        auc_sim = AUC(label_cross, presiction_similarity)
        print(auc_sim)

        Cross_Result_Log_File.write(str(epoch) + '\n')
        Cross_Result_Log_File.write('Average cross loss: ')
        Cross_Result_Log_File.write(str(ave_losss_cross) + '\n')

        Cross_Result_Log_File.write('AUC_OUT: ')
        Cross_Result_Log_File.write(str(auc) + '\n')

        Cross_Result_Log_File.write('AUC_SIM: ')
        Cross_Result_Log_File.write(str(auc_sim) + '\n')

        Cross_Result_Log_File.write('Accuracy: ')
        Cross_Result_Log_File.write(str(accuracy) + '\n')

        Cross_Result_Log_File.write('F1: ')
        Cross_Result_Log_File.write(str(F1Score) + '\n')

        Cross_Result_Log_File.flush()
        # 保存测试集的预测数据
        if auc > max(AUCList):
            saveCrossPredScore(prediction_scores, presiction_similarity, epoch)
            best_cross_model=rnn

        AccuracyList.append(accuracy)
        AUCList.append(auc)

        return prediction_scores, label_cross, auc, accuracy,best_cross_model


    def AUC(y, pred):
        fpr, tpr, thresholds = metrics.roc_curve(y, pred, pos_label=1)
        auc = metrics.auc(fpr, tpr)
        return auc


    # 保存中间结果
    def saveCrossPredScore(prediction_scores, presiction_similarity, epoch):
        epoch_pred_score_file = pred_cross_score_file + str(epoch)
        f = open(epoch_pred_score_file, 'w')
        f.write('AUC,AUC_SIM,label\n')
        for i in range(len(prediction_scores)):
            f.write(str(prediction_scores[i]) + ' ' + str(presiction_similarity[i]) + ' ' + str(label_cross[i]) + '\n')
        f.close()

    rnn = RNN()

    try:
        train(rnn)
        #torch.save(rnn, model_path)
        Cross_Result_Log_File.close()
        train_end_time = time.time()
        train_used_time = train_end_time - train_start_time
        Train_Log_File.write("\n")
        Train_Log_File.write("Start Time:" + str(time.asctime(time.localtime(train_start_time)))+"\n")
        Train_Log_File.write("End Time:" + str(time.asctime(time.localtime(train_end_time)))+"\n")
        Train_Log_File.write("Used Time:" + str(train_used_time)+"\n")
        Train_Log_File.flush()
        Train_Log_File.close()
        del rnn

    except KeyboardInterrupt:
        torch.save(rnn, model_path)
        Cross_Result_Log_File.close()
        Train_Log_File.close()
        pass




