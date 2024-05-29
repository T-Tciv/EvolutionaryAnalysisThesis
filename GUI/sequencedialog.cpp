#include "sequencedialog.h"
#include "ui_sequencedialog.h"

SequenceDialog::SequenceDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::SequenceDialog)
{
    ui->setupUi(this);
}

SequenceDialog::SequenceDialog(QWidget *parent, QString sequenceData):
    QDialog(parent),
    ui(new Ui::SequenceDialog)
{
    ui->setupUi(this);

    ui->listWidget->addItem(sequenceData);
}

SequenceDialog::~SequenceDialog()
{
    delete ui;
}
