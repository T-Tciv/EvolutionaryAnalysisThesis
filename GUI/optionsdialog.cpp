#include "optionsdialog.h"
#include "ui_optionsdialog.h"

OptionsDialog::OptionsDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::OptionsDialog)
{
    ui->setupUi(this);
    ui->programLineEdit->setText("blastn");
    ui->hitlistSizeButton->setText("1");
    ui->databaseLineEdit->setText("nt");
    ui->evalueField->setText("0.000000005");
    ui->penaltyField->setText("-2");
    ui->rewardField->setText("1");

    optionsData = new OptionsData();

}

OptionsDialog::~OptionsDialog()
{
    delete optionsData;
    delete ui;
}

void OptionsDialog::on_setOptionsButton_clicked()
{
    optionsData->setProgramName(ui->programLineEdit->text());
    optionsData->setHitlistSize(ui->hitlistSizeButton->text().toInt());
    optionsData->setDatabaseName(ui->databaseLineEdit->text());
    optionsData->setIsParsed(ui->parseCheckBox->isChecked());
    optionsData->setEvalue(ui->evalueField->text().toDouble());
    optionsData->setPenalty(ui->penaltyField->text().toDouble());
    optionsData->setReward(ui->rewardField->text().toDouble());

    emit setOptions(*optionsData);
    close();
}

