#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    blastLogic = new BlastLogic();

    // TODO: ДОБАВИТЬ ПЕРЕМЕННУЮ-СОСТОЯНИЕ, ОПРЕЛЯЮЩУЮ ТЕКУЩИЙ ВЫБОР, МЕНЯТЬ ЕЁ В СЛОТАХ
    // Ставим кнопку выбора файла в исходное состояние (блокируем)
    ui->inputRadioButton->setChecked(true);
    ui->fileSelectionButton->setDisabled(true);

    // Работа с окнами
    // Инициализируем окно-диалог для выбора файла и устаналиваем его в режим модального (т.е., пока оно открыто, больше ничего нажать нельзя)
    fileDialog = new FileSelectionDialog(this);
    fileDialog->setModal(true);
    //Получаем имя файла из этого окна по сигналу setFile
    connect(fileDialog, &FileSelectionDialog::setFile, this, &MainWindow::getFile);

    // Окно настроек:
    optionsDialog = new OptionsDialog(this);
    optionsDialog->setModal(true);
    connect(optionsDialog, &OptionsDialog::setOptions, this, &MainWindow::getOptions);

    // Окно фильтрации Entrez
    entrezDialog = new EntrezDialog(this);
    entrezDialog->setModal(true);
    connect(entrezDialog, &EntrezDialog::setEntrez, this, &MainWindow::getEntrez);

    // Поменяла кнопке поиска цвет, чтобы интерфейс был понятнее
    // TODO: вынести оформление в отдельный файл и читать из него
    ui->searchButton->setStyleSheet("background-color: #d9ead3;");

    // Поменяла метке с именем файла цвет, чтобы интерфейс был понятнее
    // TODO: вынести оформление в отдельный файл и читать из него
    ui->fileNameLabel->setStyleSheet("color: #003366;");

    // Настройка текстовых полей для флангов
    QIntValidator * intValidator = new QIntValidator();
    intValidator->setBottom(0);
    ui->leftIndentEdit->setValidator( intValidator );
    ui->rightIndentEdit->setValidator( intValidator );
    ui->leftIndentEdit->setPlaceholderText("0...");
    ui->rightIndentEdit->setPlaceholderText("0...");

    // Настройка списка последовательностей
    // TODO: ВЫНЕСТИ ОФОРМЛЕНИЕ В ОТДЕЛЬНЫЙ ФАЙЛ И ЧИТАТЬ ИЗ НЕГО
    ui->listWidget->setStyleSheet("QListWidget {padding: 10px;} QListWidget::item { border: 1px solid black; padding: 10px; margin: 10px; } QListWidget::item:selected { color: black; background-color: #d9ead6;}");
}

MainWindow::~MainWindow()
{
    delete blastLogic;
    delete entrezDialog;
    delete optionsDialog;
    delete fileDialog;
    delete ui;
}

void MainWindow::on_fileSelectionButton_clicked()
{
    fileDialog->show();
}


void MainWindow::on_listWidget_itemDoubleClicked(QListWidgetItem *item)
{
    int index = ui->listWidget->row(item);
    SequenceDialog *sequenceDialog = new SequenceDialog(this, QString::fromStdString(blastVector[index]));
    sequenceDialog->show();
}

void MainWindow::getFile(QString fileName) {
    ui->fileNameLabel->setText(fileName);
    blastLogic->setInpuFileName(fileName.toStdString());

    ifstream fastaFile (fileName.toStdString());

    std::string fastaFileContent ((std::istreambuf_iterator<char>(fastaFile)), std::istreambuf_iterator<char>());

    std::string delimiter = ">";

    size_t pos = 0;
    std::string sequenceText;
    while ((pos = fastaFileContent.find(delimiter)) != std::string::npos) {
        sequenceText = fastaFileContent.substr(0, pos);

        if (pos != 0) {
            ui->listWidget->addItem(QString::fromStdString(sequenceText));
        }

        fastaFileContent.erase(0, pos + delimiter.length());
    }

    ui->listWidget->addItem(QString::fromStdString(fastaFileContent));
}

void MainWindow::getOptions(OptionsData optionsSet) {
    blastLogic->setBlastSettings(optionsSet.getProgramName().toStdString(), optionsSet.getHitlistSize(), optionsSet.getDatabaseName().toStdString(), optionsSet.getIsParsed());
    blastLogic->setEvalue(optionsSet.getEvalue());
    blastLogic->setReward(optionsSet.getReward());
    blastLogic->setPenalty(optionsSet.getPenalty());
}

void MainWindow::getEntrez(QString entrezQuery) {
    blastLogic->setEntrezQuery(entrezQuery.toStdString());
}


void MainWindow::on_pushButton_clicked()
{
    optionsDialog->show();
}


void MainWindow::on_inputRadioButton_toggled(bool checked)
{
    QTextEdit *inputFiled = ui->inputField;
    inputFiled->setDisabled(!checked);
}


void MainWindow::on_fileRadioButton_toggled(bool checked)
{
    QPushButton *fileButton = ui->fileSelectionButton;
    fileButton->setDisabled(!checked);
}


void MainWindow::on_filterButton_clicked()
{
    entrezDialog->show();
}


void MainWindow::on_searchButton_clicked()
{
    blastLogic->setIndents(ui->leftIndentEdit->text().toInt(), ui->rightIndentEdit->text().toInt());
    blastVector = blastLogic->makeBlast();
}

