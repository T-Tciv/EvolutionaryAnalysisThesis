#ifndef SEQUENCEDIALOG_H
#define SEQUENCEDIALOG_H

#include <QDialog>

namespace Ui {
class SequenceDialog;
}

class SequenceDialog : public QDialog
{
    Q_OBJECT

public:
    explicit SequenceDialog(QWidget *parent = nullptr);
    SequenceDialog(QWidget *parent = nullptr, QString sequenceData = "");
    ~SequenceDialog();

private:
    Ui::SequenceDialog *ui;
};

#endif // SEQUENCEDIALOG_H
