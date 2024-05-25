function [accuracy] = SVM_Trainer(O, Y)
    % Concatenate data and labels
    data = [O; Y];
    labels = [ones(size(O, 1), 1); zeros(size(Y, 1), 1)];

    % Train SVM classifier
    svm_model = fitcsvm(data, labels);

    % Predict labels for training data
    predicted_labels_train = predict(svm_model, data);

    % Evaluate training accuracy
    train_accuracy = mean(predicted_labels_train == labels) * 100;
    fprintf('Training Accuracy: %.2f%%\n', train_accuracy);

    % Cross-validation accuracy (using training data)
    accuracy = train_accuracy;
end
