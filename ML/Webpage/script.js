// const PRECISE = 0.9;
const SENSITIVE = 0.5;

function openTab(tabName, position) {
    let i;
    let x = document.getElementsByClassName("tab");
    for (i = 0; i < x.length; i++) {
        x[i].style.display = "none";
    }
    document.getElementById(tabName).style.display = "block";
    // document.querySelector('.bar-button.button1.selected').className = "bar-button button1";
    // document.getElementById(`${tabName}Tab`).className = "bar-button button1 selected";

    document.getElementById('navTool').style.left = ((100 / document.getElementById('menuItems').childElementCount) * position) + '%';
}


function simpleSeq(x) {
    return {
        name: '',
        seq: x
    };
}

function loadExample(tab) {
    switch (tab) {
        case 'miRNA':
            document.querySelector(`#${tab} #smallRNA`).value = "TCCGAGCCTGGGTCTCCCTC\nAAAGTGCTTCCCTTTGGACT";
            document.querySelector(`#${tab} #target`).value = "TTCAGGAGAAGCTGAGAGAGACCCAGGAGTATAACCGAATTCAGAAGGAG\nCCTGAGGAGACCAAGCTGGGCAAGAGGCAGTGCGACGGCAAGAATGCGCT";
            return;
        case 'tRNA':
            document.querySelector(`#${tab} #smallRNA`).value = "TCCCTGGTGGTCTAGTGGTT\nAACCAGGCGGAAACACCAAA";
            document.querySelector(`#${tab} #target`).value = "GACCACCACCTCAGCTCTGCGGACCTTTTGGGCCTCGGCCACTTCCTCCA\nTCCACCAGGAAGACGGGACCTGCCTCTCCACCCTCGGGGATTTTTACCTG";
            return;
    }
}

function getSeq(tab) {
    let error = false;
    let warning = [];

    const validate = (array, field = 'small RNA', limit = 20) => {
        array.map((item, index) => {
            item.toUpperCase().replace(/\s/g, '').replace(/U/g, 'T');
            if (!item.length) {
                error = true;
                warning.push(`Error: ${field} on line: ${index + 1} is empty`);
            } else if (item.length > limit) {
                warning.push(`Warning: ${field} on line: ${index + 1} is longer than ${limit} bp. Trimming.`);
                item = item.slice(0, limit);
            } else if (item.length < limit) {
                warning.push(`Warning: ${field} on line: ${index + 1} is shorter than ${limit} bp. Too short sequences may lead to unexpected behaviour.`);
            }

            if (!/^([ACGT])*$/.test(item)) {
                error = true;
                warning.push(`Error: ${field} on line: ${index + 1} contains not allowed bases. Allowed bases are: "A", "C", "G", "T", "U"`);
            }

            return item;
        });

        return array;
    };
    const smallRNA = validate(document.querySelector(`#${tab} #smallRNA`).value?.trim()?.split(/\r?\n/));
    const target = validate(document.querySelector(`#${tab} #target`).value?.trim()?.split(/\r?\n/), 'target', 50);

    if (smallRNA.length !== target.length) {
        error = true;
        warning.push(`Error: The number of miRNAs and mRNAs must be equal.`);
    }

    return {
        smallRNA: smallRNA,
        target: target,
        error: error,
        warning: warning,
    };
}

async function makePrediction(tab) {
    const seq = getSeq(tab);
    document.querySelector(`#${tab} #result`).innerHTML = 'Processing...';
    document.querySelector(`#${tab} #warning`).innerHTML = seq.warning.join('<br />');

    if (seq.error) {
        return;
    }

    const resultColor = (x) => {
        switch (true) {
            // case x > PRECISE:
            //     return 'blue';
            case x > SENSITIVE:
                return 'blue';
            default:
                return '';
        }
    }
    let results = [];
    let modelUrl = getModelUrl(tab);
    document.querySelector(`#${tab} #result`).innerHTML += '<br>Loading model...';
    let model = await tf.loadLayersModel(modelUrl);
    document.querySelector(`#${tab} #result`).innerHTML += '<br>Model loaded<br>Processing data...';
    for (var i = 0; i < seq.smallRNA.length; i++) {
        const smallRNA = seq.smallRNA[i];
        const target = seq.target[i];
        let input = "";
        for (var j = 0; j < 1000; j++) {
            input[j] = new Array(20);
        }
        binding = ["AT", "TA", "CG", "GC"];
        for (var k = 0; k < 50; k++) {
            for (var l = 0; l < 20; l++) {
                char = smallRNA.charAt(l) + target.charAt(k);
                if (char.length !== 2) {
                    input += '0';
                } else if (binding.includes(char)) {
                    input += '1';
                } else {
                    input += '0';
                }
            }
        }
        y1 = tf.tensor1d(input.split(''), 'int32');
        y2 = y1.reshape([1, 50, 20, 1]);

        let result = parseFloat(model.predict(y2)?.asScalar()?.dataSync())?.toFixed(4);
        results.push("<tr><td>" + smallRNA + "</td><td>" + target + "</td><td style='color:" + resultColor(result) + ";'>" + result + "</td></tr>");
    }
    document.querySelector(`#${tab} #result`).innerHTML = '<table><thead><th>miRNA sequence</th><th>mRNA sequence</th><th>score</th></thead>' +
        results.join('') +
        '</table>'
}

function getModelUrl(tab) {
    switch (tab) {
        case 'miRNA':
            return "https://raw.githubusercontent.com/evaklimentova/smallRNA_binding/main/Webpage/tfjs_miRNA_model/model.json";
        case 'tRNA':
            return "https://raw.githubusercontent.com/evaklimentova/smallRNA_binding/main/Webpage/tfjs_tRNA_model/model.json";
    }
}
