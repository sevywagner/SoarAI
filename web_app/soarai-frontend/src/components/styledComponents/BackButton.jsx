import { styled } from "styled-components";

const BackButton = styled.button`
    border: 1px solid white;
    background-color: black;
    color: white;
    font-size: 2.5vh;
    margin: 2vh;
    width: 100%;
    padding: 2vh;
    border-radius: 3vh;

    &:hover {
        animation: hover .2s ease-in-out;
        background-color: white;
        color: black;
    }

    @keyframes hover {
        0% {
            border: 1px solid white;
            background-color: black;
            color: white;
        }

        100% {
            background-color: white;
            color: black;
        }
    }
`

export default BackButton;