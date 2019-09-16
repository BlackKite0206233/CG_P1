/********************************************************************************
** Form generated from reading UI file 'prompt.ui'
**
** Created by: Qt User Interface Compiler version 5.13.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_PROMPT_H
#define UI_PROMPT_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QDialog>
#include <QtWidgets/QLabel>
#include <QtWidgets/QPushButton>

QT_BEGIN_NAMESPACE

class Ui_Dialog
{
public:
    QPushButton *GLButton;
    QPushButton *DirectXButton;
    QLabel *label;

    void setupUi(QDialog *Dialog)
    {
        if (Dialog->objectName().isEmpty())
            Dialog->setObjectName(QString::fromUtf8("Dialog"));
        Dialog->setWindowModality(Qt::ApplicationModal);
        Dialog->resize(401, 228);
        Dialog->setContextMenuPolicy(Qt::NoContextMenu);
        Dialog->setSizeGripEnabled(true);
        Dialog->setModal(true);
        GLButton = new QPushButton(Dialog);
        GLButton->setObjectName(QString::fromUtf8("GLButton"));
        GLButton->setGeometry(QRect(20, 80, 361, 51));
        DirectXButton = new QPushButton(Dialog);
        DirectXButton->setObjectName(QString::fromUtf8("DirectXButton"));
        DirectXButton->setGeometry(QRect(20, 140, 361, 51));
        label = new QLabel(Dialog);
        label->setObjectName(QString::fromUtf8("label"));
        label->setGeometry(QRect(10, 20, 371, 20));
        label->setLayoutDirection(Qt::RightToLeft);
        label->setAlignment(Qt::AlignCenter);

        retranslateUi(Dialog);

        QMetaObject::connectSlotsByName(Dialog);
    } // setupUi

    void retranslateUi(QDialog *Dialog)
    {
        Dialog->setWindowTitle(QCoreApplication::translate("Dialog", "Dialog", nullptr));
        GLButton->setText(QCoreApplication::translate("Dialog", "OpenGL", nullptr));
        DirectXButton->setText(QCoreApplication::translate("Dialog", "DirectX", nullptr));
        label->setText(QCoreApplication::translate("Dialog", "Choose Render", nullptr));
    } // retranslateUi

};

namespace Ui {
    class Dialog: public Ui_Dialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_PROMPT_H
